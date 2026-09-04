// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

using namespace Rcpp;

#define PACK_DENSITY 4

/* 3 is 11 in binary, we need a 2 bit mask for each of the 4 positions */
#define MASK0 3	  /* 3 << 2 * 0 */
#define MASK1 12  /* 3 << 2 * 1 */
#define MASK2 48  /* 3 << 2 * 2 */
#define MASK3 192 /* 3 << 2 * 3 */


// ---------------------------------------------------------------------------
// packedancestrymap (.geno) layout
//
// Both layouts described in DReichLab/AdmixTools transpose_info/transpose_format.txt
// start with a whitespace-delimited text header
//
//     <GENO|TGENO> <nind> <nsnp> <indhash> <snphash>
//
// followed by a byte-aligned matrix of 2-bit genotype codes (0, 1, 2 as in
// the EIGENSTRAT genotype value, 3 = missing), four per byte, first code
// in the two high bits.
//
//   GENO   (packed ancestrymap, the original layout): the header occupies
//          max(48, ceil(nind/4)) bytes; then one row per SNP of the same
//          size, holding the genotypes of all individuals.
//   TGENO  (transpose packed, what recent AADR releases ship): the header
//          occupies exactly 48 bytes; then one row per INDIVIDUAL of
//          ceil(nsnp/4) bytes, holding that individual's genotypes. The
//          reference writer does not pad these rows, but plink2's
//          `--export eigt` pads them to 48 bytes, so both are accepted.
//
// The previous implementation inferred the row size as filesize/(nsnp+1)
// and always read SNP-major rows, which silently returns chance-level
// genotypes for a TGENO file and, when the inferred size rounds down,
// writes past the end of its per-row buffers (uqrmaie1/admixtools#172).
// The readers below take the layout from the header and check that the
// file size matches it before any genotype byte is interpreted.
// ---------------------------------------------------------------------------

struct GenoLayout {
  bool transposed;        // false: GENO (row per SNP); true: TGENO (row per individual)
  long long nind, nsnp;   // as declared in the header
  long long headerbytes;  // offset of the first matrix byte
  long long rowbytes;     // bytes per matrix row
};

static GenoLayout read_geno_layout(std::ifstream& in, const std::string& path, int nsnp, int nind,
                                   const IntegerVector& indvec) {

  if(nsnp < 0 || nind < 0) {
    stop("Negative number of SNPs (%d) or individuals (%d) requested for %s", nsnp, nind, path);
  }
  if(indvec.size() != nind) {
    stop("indvec has %d entries but %d individuals were expected for %s "
         "(does the .ind file match the genotype file?)", (int)indvec.size(), nind, path);
  }
  in.seekg(0, std::ifstream::end);
  long long len = (long long)in.tellg();
  if(len < 48) {
    stop("Genotype file %s has only %d bytes, too small for a packedancestrymap header", path, len);
  }
  char hdr[49] = {0};
  in.seekg(0, std::ifstream::beg);
  if(!in.read(hdr, 48)) {
    stop("Could not read the header of genotype file %s", path);
  }
  std::istringstream ss(hdr);
  std::string tag;
  long long hind = -1, hsnp = -1;
  ss >> tag >> hind >> hsnp;

  GenoLayout g;
  if(tag == "GENO") {
    g.transposed = false;
  } else if(tag == "TGENO") {
    g.transposed = true;
  } else {
    stop("Genotype file %s does not start with a 'GENO' (packedancestrymap) or 'TGENO' "
         "(transposed packedancestrymap) header", path);
  }
  if(hind != nind || hsnp != nsnp) {
    stop("Header of %s declares %d individuals and %d SNPs, but the caller expected %d individuals "
         "and %d SNPs (usually the row counts of the .ind and .snp files)", path, hind, hsnp, nind, nsnp);
  }
  g.nind = hind;
  g.nsnp = hsnp;
  long long nrows;
  if(g.transposed) {
    g.headerbytes = 48;
    g.rowbytes = (g.nsnp + 3) / 4;
    nrows = g.nind;
    // plink2 pads short TGENO rows to the 48-byte header size.
    if(g.rowbytes < 48 && len != g.headerbytes + g.rowbytes * nrows && len == g.headerbytes + 48 * nrows) {
      g.rowbytes = 48;
    }
  } else {
    g.rowbytes = std::max((long long)48, (g.nind + 3) / 4);
    g.headerbytes = g.rowbytes;
    nrows = g.nsnp;
  }
  long long expected = g.headerbytes + g.rowbytes * nrows;
  if(len != expected) {
    stop("Genotype file %s has %d bytes, but its %s header (%d individuals, %d SNPs) implies %d bytes",
         path, len, tag, g.nind, g.nsnp, expected);
  }
  return g;
}

// 2-bit genotype code of element k (0-based) in a packed row.
static inline int unpack2(const unsigned char* row, long long k) {
  return (row[k >> 2] >> (6 - 2 * (k & 3))) & 3;
}


// [[Rcpp::export]]
IntegerMatrix cpp_read_packedancestrymap(String genofile, int nsnp, int nind, IntegerVector indvec,
                                         int first, int last, bool transpose = false, bool verbose = true) {

  // Returns an IntegerMatrix of genotype dosages in {0, 1, 2, NA_INTEGER},
  // matching cpp_read_plink. Genotypes are integer-valued; an integer
  // return halves the per-block buffer and lets the block feed
  // cpp_gmat_to_aftable (arma::imat) without a double->int coercion copy.
  //
  // `transpose` only chooses the shape of the returned matrix
  // (individuals x SNPs instead of SNPs x individuals). The on-disk layout
  // (GENO or TGENO) is taken from the file header, see read_geno_layout.
  int val;
  long long bytespersnp;
  int readsnps = last - first;

  std::ifstream in(genofile.get_cstring(), std::ios::in | std::ios::binary);

  if(!in) {
    stop("Genotype file not found!");
  }
  GenoLayout g = read_geno_layout(in, genofile.get_cstring(), nsnp, nind, indvec);
  if(first < 0 || last > nsnp || readsnps < 0) {
    stop("SNP range [%d, %d) is outside the %d SNPs in the genotype file", first, last, nsnp);
  }

  int nindused = 0;
  for(int i = 0 ; i < nind; i++) {
    if(indvec[i] == 1) nindused++;
  }

  IntegerMatrix geno(transpose?nindused:readsnps, transpose?readsnps:nindused);
  std::fill(geno.begin(), geno.end(), NA_INTEGER);

  if(g.transposed) {
    // TGENO: one row per individual. For each kept individual read only
    // the bytes spanning SNPs [first, last) of its row.
    if(readsnps > 0) {
      long long firstbyte = first / PACK_DENSITY;
      long long nbytes = (last - 1) / PACK_DENSITY - firstbyte + 1;
      std::vector<unsigned char> row(nbytes);
      int c = 0;
      for(int i = 0; i < nind; i++) {
        if(!indvec[i]) continue;
        if(verbose && c % 1000 == 0) Rcout << "\r" << c/1000 << "k samples read...";
        in.seekg(g.headerbytes + (long long)i * g.rowbytes + firstbyte, std::ifstream::beg);
        in.read((char*)row.data(), nbytes);
        for(int j = 0; j < readsnps; j++) {
          val = unpack2(row.data(), (long long)first + j - PACK_DENSITY * firstbyte);
          if(val != 3) {
            if(transpose) geno(c, j) = val;
            else geno(j, c) = val;
          }
        }
        c++;
      }
    }
    if(verbose) Rcout << std::endl;
    in.close();
    return geno;
  }

  // GENO: one row per SNP of g.rowbytes bytes, header of the same size.
  bytespersnp = g.rowbytes;

  std::vector<int> blockused(bytespersnp, 0);
  for(int i = 0 ; i < nind; i++) {
    if(indvec[i] == 1) {
      blockused[i/PACK_DENSITY] = 1;
    }
  }

  in.seekg(g.headerbytes + (long long)first * bytespersnp, std::ifstream::beg);
  std::vector<char> tmp(bytespersnp + 1, 0);
  char tmpi;

  // Allocate more than the sample size since data must take up whole bytes
  std::vector<char> tmp2(bytespersnp * PACK_DENSITY + 1, 0);

  int k;
  for(int j = 0 ; j < readsnps; j++) {
    //for(unsigned int j = 0 ; j < 3; j++) {
    if(verbose && j % 1000 == 0) Rcout << "\r" << j/1000 << "k SNPs read...";

    // read raw genotypes
    in.read(tmp.data(), sizeof(char) * bytespersnp);

    for(int l = 0; l < bytespersnp; l++) {
      if(!blockused[l]) continue;

      tmpi = tmp[l];
      k = PACK_DENSITY * l;

      /* geno is interpreted as a char, however a1 and a2 are bits for allele 1 and
       * allele 2. The final genotype is the sum of the alleles, except for 11
       * which denotes missing.
       */
      tmp2[k] = (tmpi & MASK3) >> 6;
      tmp2[k+1] = (tmpi & MASK2) >> 4;
      tmp2[k+2] = (tmpi & MASK1) >> 2;
      tmp2[k+3] = (tmpi & MASK0);
    }

    int c = 0;
    if(!transpose) {
      for(int i = 0; i < nind; i++) {
        if(!indvec[i]) continue;
        val = tmp2[i];
        if(val != 3) geno(j, c) = val;
        c++;
      }
    } else {
      for(int i = 0; i < nind; i++) {
        if(!indvec[i]) continue;
        val = tmp2[i];
        if(val != 3) geno(c, j) = val;
        c++;
      }
    }
  }
  if(verbose) Rcout << std::endl;

  in.close();

  return geno;
}


// [[Rcpp::export]]
NumericVector cpp_packedancestrymap_ploidy(String genofile, int nsnp, int nind, IntegerVector indvec, int ntest = 1000) {

  int val, k;
  long long bytespersnp;

  if(ntest == NA_INTEGER || ntest < 0) {
    stop("ntest must be a non-negative integer");
  }
  ntest = std::min(ntest, nsnp);
  std::ifstream in(genofile.get_cstring(), std::ios::in | std::ios::binary);
  if(!in) {
    stop("Genotype file not found!");
  }
  GenoLayout g = read_geno_layout(in, genofile.get_cstring(), nsnp, nind, indvec);

  NumericVector ploidy(nind);
  ploidy.fill(1.0);

  if(g.transposed) {
    // TGENO: scan the first ntest SNPs of each used individual's row.
    long long nbytes = (ntest + PACK_DENSITY - 1) / PACK_DENSITY;
    std::vector<unsigned char> row(nbytes);
    for(int i = 0; i < nind; i++) {
      if(indvec[i] == -1) continue;
      in.seekg(g.headerbytes + (long long)i * g.rowbytes, std::ifstream::beg);
      in.read((char*)row.data(), nbytes);
      for(int j = 0; j < ntest; j++) {
        if(unpack2(row.data(), j) == 1) {
          ploidy(i) = 2.0;
          break;
        }
      }
    }
    in.close();
    return ploidy;
  }

  bytespersnp = g.rowbytes;
  std::vector<int> blockused(bytespersnp, 0);
  for(int i = 0 ; i < nind; i++) {
    if(indvec[i] != -1) {
      blockused[i/PACK_DENSITY] = 1;
    }
  }
  std::vector<unsigned char> tmp(bytespersnp, 0);
  unsigned char tmpi;
  std::vector<unsigned char> tmp2(bytespersnp * PACK_DENSITY, 0);

  in.seekg(g.headerbytes, std::ifstream::beg);
  for(int j = 0 ; j < ntest; j++) {
    in.read((char*)tmp.data(), sizeof(char) * bytespersnp);

    for(int l = 0 ; l < bytespersnp; ++l) {
      if(!blockused[l]) continue;
      tmpi = tmp[l];
      k = PACK_DENSITY * l;
      tmp2[k] = (tmpi & MASK3) >> 6;
      tmp2[k+1] = (tmpi & MASK2) >> 4;
      tmp2[k+2] = (tmpi & MASK1) >> 2;
      tmp2[k+3] = (tmpi & MASK0);
    }
    for(int i = 0; i < nind; i++) {
      if(indvec[i] == -1) continue;
      val = (double)tmp2[i];
      if(val == 1) ploidy(i) = 2.0;
    }
  }

  in.close();

  return ploidy;
}


// [[Rcpp::export]]
List cpp_packedancestrymap_to_afs(String genofile, int nsnp, int nind, IntegerVector indvec,
                                      int first, int last, IntegerVector ploidy,
                                      bool transpose, bool verbose) {
  //same arguments as cpp_read_packedancestrymap, except indvec assigns populations. indiv not used: -1

  int val, k, pop;

  long long bytespersnp;
  int readsnps = std::min(nsnp, last?last:nsnp) - first;

  std::ifstream in(genofile.get_cstring(), std::ios::in | std::ios::binary);

  if(!in) {
    stop("Genotype file not found!");
  }
  GenoLayout g = read_geno_layout(in, genofile.get_cstring(), nsnp, nind, indvec);
  if(first < 0 || readsnps < 0) {
    stop("SNP range starting at %d is outside the %d SNPs in the genotype file", first, nsnp);
  }
  if(ploidy.size() != nind) {
    stop("ploidy has %d entries but %d individuals were expected", (int)ploidy.size(), nind);
  }

  int numpop = max(indvec)+1;

  NumericMatrix afs(transpose?numpop:readsnps, transpose?readsnps:numpop);
  std::fill(afs.begin(), afs.end(), NA_REAL);
  NumericMatrix counts(transpose?numpop:readsnps, transpose?readsnps:numpop);
  std::fill(counts.begin(), counts.end(), 0.0);

  if(g.transposed) {
    // TGENO: one row per individual. Accumulate allele sums directly into
    // `afs` and observed alleles into `counts`, then normalise in place.
    // Contributions to each cell are added in ascending individual order,
    // as in the GENO branch, so both layouts give identical numbers.
    if(readsnps > 0) {
      std::fill(afs.begin(), afs.end(), 0.0);
      long long firstbyte = first / PACK_DENSITY;
      long long nbytes = ((long long)first + readsnps - 1) / PACK_DENSITY - firstbyte + 1;
      std::vector<unsigned char> row(nbytes);
      int c = 0;
      for(int i = 0; i < nind; i++) {
        if(indvec[i] == -1) continue;
        if(verbose && c % 1000 == 0) Rcout << "\r" << c/1000 << "k samples read...";
        c++;
        pop = indvec(i);
        in.seekg(g.headerbytes + (long long)i * g.rowbytes + firstbyte, std::ifstream::beg);
        in.read((char*)row.data(), nbytes);
        for(int j = 0; j < readsnps; j++) {
          val = unpack2(row.data(), (long long)first + j - PACK_DENSITY * firstbyte);
          if(val != 3) {
            if(!transpose) {
              afs(j, pop) += val / (3.0-ploidy(i));
              counts(j, pop) += ploidy(i);
            } else {
              afs(pop, j) += val / (3.0-ploidy(i));
              counts(pop, j) += ploidy(i);
            }
          }
        }
      }
      for(int j = 0; j < readsnps; j++) {
        for(int p = 0; p < numpop; p++) {
          if(!transpose) {
            if(counts(j, p) > 0) afs(j, p) = afs(j, p) / counts(j, p);
            else afs(j, p) = NA_REAL;
          } else {
            if(counts(p, j) > 0) afs(p, j) = afs(p, j) / counts(p, j);
            else afs(p, j) = NA_REAL;
          }
        }
      }
    }
    if(verbose) Rcout << std::endl;
    in.close();
    return Rcpp::List::create(_["afs"] = afs, _["counts"] = counts);
  }

  bytespersnp = g.rowbytes;

  std::vector<int> blockused(bytespersnp, 0);
  for(int i = 0 ; i < nind; i++) {
    if(indvec[i] != -1) {
      blockused[i/PACK_DENSITY] = 1;
    }
  }

  //set ploidy
  std::vector<unsigned char> tmp(bytespersnp, 0);
  unsigned char tmpi;
  // Allocate more than the sample size since data must take up whole bytes
  std::vector<unsigned char> tmp2(bytespersnp * PACK_DENSITY, 0);

  in.seekg(g.headerbytes + (long long)first * bytespersnp, std::ifstream::beg);

  NumericVector altalleles(numpop);
  NumericVector observedalleles(numpop);

  for(int j = 0 ; j < readsnps; j++) {
    if(verbose && j % 1000 == 0) Rcout << "\r" << j/1000 << "k SNPs read...";

    // read raw genotypes
    in.read((char*)tmp.data(), sizeof(char) * bytespersnp);

    for(int l = 0 ; l < bytespersnp; ++l) {
      if(!blockused[l]) continue;

      tmpi = tmp[l];
      k = PACK_DENSITY * l;

      /* geno is interpreted as a char, however a1 and a2 are bits for allele 1 and
       * allele 2. The final genotype is the sum of the alleles, except for 11
       * which denotes missing.
       */
      tmp2[k] = (tmpi & MASK3) >> 6;
      tmp2[k+1] = (tmpi & MASK2) >> 4;
      tmp2[k+2] = (tmpi & MASK1) >> 2;
      tmp2[k+3] = (tmpi & MASK0);
    }

    altalleles.fill(0);
    observedalleles.fill(0);
    for(int i = 0; i < nind; i++) {
      if(indvec[i] == -1) continue;
      val = (double)tmp2[i];
      if(val != 3) {
        pop = indvec(i);
        altalleles(pop) += val / (3.0-ploidy(i));
        observedalleles(pop) += ploidy(i);
      }
    }

    if(!transpose) {
      for(int i = 0; i < numpop; i++) {
        if(observedalleles(i) > 0) {
          afs(j, i) = altalleles(i)/observedalleles(i);
          counts(j, i) = observedalleles(i);
        }
      }
    } else {
      for(int i = 0; i < numpop; i++) {
        if(observedalleles(i) > 0) {
          afs(i, j) = altalleles(i)/observedalleles(i);
          counts(i, j) = observedalleles(i);
        }
      }
    }
  }
  if(verbose) Rcout << std::endl;

  in.close();

  return Rcpp::List::create(_["afs"] = afs, _["counts"] = counts);
}


// [[Rcpp::export]]
IntegerMatrix cpp_read_eigenstrat(String genofile, int nsnp, int nind, IntegerVector indvec,
                                  int first, int last, bool transpose = false, bool verbose = true) {

  // Returns an IntegerMatrix of genotype dosages in {0, 1, 2, NA_INTEGER},
  // matching cpp_read_plink. See that function for the rationale; the
  // missing sentinel in EIGENSTRAT is 9 (vs 3 in BED/PACKEDANCESTRYMAP).
  int val;
  long long len, bytespersnp;
  int readsnps = last - first;

  std::ifstream in(genofile.get_cstring(), std::ios::in | std::ios::binary);

  if(!in) {
    stop("Genotype file not found!");
  }
  in.seekg(0, std::ifstream::end);
  // file size in bytes
  len = (long long)in.tellg();
  bytespersnp = len/nsnp;

  int nindused = 0;
  for(int i = 0 ; i < nind; i++) {
    if(indvec[i] == 1) {
      nindused++;
    }
  }

  IntegerMatrix geno(transpose?nindused:readsnps, transpose?readsnps:nindused);
  std::fill(geno.begin(), geno.end(), NA_INTEGER);
  in.seekg(first*bytespersnp, std::ifstream::beg);
  std::string tmp;
  for(int j = 0 ; j < readsnps; j++) {
    int c = 0;
    if(verbose && j % 1000 == 0) Rcout << "\r" << j/1000 << "k SNPs read...";
    std::getline(in, tmp);
    if(!transpose) {
      for(int i = 0; i < nind; i++) {
        if(!indvec[i]) continue;
        val = tmp[i]-'0';
        if(val != 9) geno(j, c) = val;
        c++;
      }
    } else {
      for(int i = 0; i < nind; i++) {
        if(!indvec[i]) continue;
        val = tmp[i]-'0';
        if(val != 9) geno(c, j) = val;
        c++;
      }
    }
  }
  if(verbose) Rcout << std::endl;
  in.close();
  return geno;
}

