// Oracle harness for the Rust port -- NOT part of CoCo.
//
// Several CountProfile functions are unreachable from the command line: some are
// private, one is a file-static helper, and two (checkForSpuriousTransitionDrops
// and checkForSpuriousTransitionDropsWithWindow) have no caller at all in the
// released code. End-to-end diffing therefore cannot pin them. This binary calls
// them directly on synthetic inputs and prints their results as text, so the Rust
// tests can compare function by function rather than only file by file.
//
// CountProfile.cpp is #included rather than linked so that the file-static
// calcNeighborhoodTolerance is reachable, and `private` is redefined so the
// profile array can be set without going through fill(). The whole standard
// library and all of CoCo's other headers are pulled in *before* that redefinition
// -- otherwise libstdc++ headers get parsed with `private` meaning `public`, and
// classes that declare the same member twice under different access (std::
// stringbuf does) fail to compile.

#include <bits/stdc++.h>
#include <gatb/gatb_core.hpp>

// CoCo's own headers are pulled in *by* CountProfile.cpp, under the redefinition,
// so that KmerTranslator's mask arrays and CountProfile's profile array are both
// reachable. Including them beforehand would fix their access as private.
#define private public
#include "../../src/CountProfile.cpp"
#undef private

// kmer.h declares minIndex and packedKmer2String but not revComplement, which is
// nonetheless an external symbol in kmer.cpp.
packedKmerType revComplement(packedKmerType kmer, unsigned short kmerSize);

/* A LookupTableBase backed by std::map, so a test can state a table exactly
   instead of going through DSK. */
class MapTable : public LookupTableBase {
public:
  std::map<packedKmerType, unsigned int> m;
  unsigned int getCount(const packedKmerType kmer) const {
    auto it = m.find(kmer);
    return it == m.end() ? 0u : it->second;
  }
  bool increaseCount(packedKmerType kmer) {
    auto it = m.find(kmer);
    if (it == m.end() || it->second == 0) return false;
    it->second += 1; return true;
  }
  bool decreaseCount(packedKmerType kmer) {
    auto it = m.find(kmer);
    if (it == m.end() || it->second == 0) return false;
    it->second -= 1; return true;
  }
  void iterateOverAll(FILE *fp) const {}
  void printSize() {}
};

/* Count every spaced k-mer of `ref` with weight `mult`, both strands via minIndex. */
static void buildTableFromReference(MapTable &t, const KmerTranslator &tr,
                                    const std::string &ref, unsigned int mult) {
  unsigned short span = tr.getSpan();
  if (ref.size() < span) return;
  spacedKmerType kmer = 0, nStore = 0;
  for (size_t idx = 0; idx < ref.size(); idx++) {
    char code = (char) res2int[(int) ref[idx]];
    if (code != -1) { kmer = (kmer << 2) | code; nStore = nStore << 1; }
    else { kmer = kmer << 2; nStore = nStore << 1 | 1; }
    if ((int) idx >= span - 1) {
      if ((nStore & tr._spaced_mask) != 0) continue;
      t.m[tr.kmer2minPackedKmer(kmer)] += mult;
    }
  }
}

static std::vector<long> parseInts(const std::string &s) {
  std::vector<long> v;
  std::istringstream is(s);
  long x;
  while (is >> x) v.push_back(x);
  return v;
}

/* Install a profile directly: value -1 marks an invalid position. */
static void setProfile(CountProfile &cp, const std::vector<long> &vals,
                       const KmerTranslator &tr) {
  delete[] cp.profile;
  cp.profile = new CountProfileEntry[vals.size()];
  for (size_t i = 0; i < vals.size(); i++) {
    cp.profile[i].valid = vals[i] < 0 ? 0 : 1;
    cp.profile[i].count = vals[i] < 0 ? 0 : (uint32_t) vals[i];
    cp.profile[i].kmer = 0;
  }
  cp.profile_length = vals.size();
  cp.profile_length_alloc = vals.size();
  cp.translator = &tr;
}

int main(int argc, char **argv) {
  initialize();
  if (argc < 2) { fprintf(stderr, "usage: oracle <mode> [args]\n"); return 2; }
  std::string mode = argv[1];
  std::string line;

  if (mode == "revcomp" || mode == "minindex" || mode == "packed2string") {
    unsigned short k = atoi(argv[2]);
    while (std::getline(std::cin, line)) {
      if (line.empty()) continue;
      uint64_t v = strtoull(line.c_str(), NULL, 16);
      if (mode == "packed2string") {
        char *s = packedKmer2String(v, k);
        printf("%s\n", s);
        free(s);
      } else {
        uint64_t r = mode == "revcomp" ? revComplement(v, k) : minIndex(v, k);
        printf("%016llx\n", (unsigned long long) r);
      }
    }
    return 0;
  }

  if (mode == "translator") {
    KmerTranslator tr(argv[2]);
    unsigned int li, lo;
    tr.getBestSplit(li, lo);
    printf("span\t%u\n", tr.getSpan());
    printf("weight\t%u\n", tr.getWeight());
    printf("longestBlock\t%u\n", tr.getLongestBlock());
    printf("logIndexSize\t%u\n", li);
    printf("logOffsetSize\t%u\n", lo);
    printf("maskArray");
    for (unsigned i = 0; i < tr.getWeight(); i++) printf("\t%u", tr._mask_array[i]);
    printf("\n");
    printf("inverseMaskArray");
    for (unsigned i = 0; i < tr.getSpan(); i++) printf("\t%u", tr._inverse_mask_array[i]);
    printf("\n");
    return 0;
  }

  if (mode == "avgqual") {
    while (std::getline(std::cin, line)) {
      if (line.empty()) continue;
      size_t sp = line.rfind(' ');
      std::string q = line.substr(0, sp);
      unsigned int pos = atoi(line.c_str() + sp + 1);
      printf("%d\n", (int) (unsigned char) getAvgQual(q, pos));
    }
    return 0;
  }

  if (mode == "maximize") {
    KmerTranslator tr(argv[2]);
    MapTable table;
    CountProfile cp(&tr, &table);
    while (std::getline(std::cin, line)) {
      if (line.empty()) continue;
      setProfile(cp, parseInts(line), tr);
      uint32_t *mx = cp.maximize();
      size_t n = cp.profile_length + tr.getSpan() - 1;
      for (size_t i = 0; i < n; i++) printf(i ? " %u" : "%u", mx[i]);
      printf("\n");
      delete[] mx;
    }
    return 0;
  }

  if (mode == "tolerance") {
    KmerTranslator tr(argv[2]);
    double threshold = atof(argv[3]);
    unsigned int pseudo = atoi(argv[4]), lower = atoi(argv[5]);
    while (std::getline(std::cin, line)) {
      if (line.empty()) continue;
      std::vector<long> v = parseInts(line);
      std::vector<uint32_t> mp(v.begin(), v.end());
      std::vector<uint32_t> tol(mp.size());
      calcNeighborhoodTolerance(tol.data(), mp.data(), mp.size(), tr.getSpan(),
                                threshold, pseudo, lower);
      for (size_t i = 0; i < tol.size(); i++) printf(i ? " %u" : "%u", tol[i]);
      printf("\n");
    }
    return 0;
  }

  if (mode == "quantile") {
    double q = atof(argv[2]);
    KmerTranslator tr("11110111111011011101010111011011111101111");
    MapTable table;
    CountProfile cp(&tr, &table);
    while (std::getline(std::cin, line)) {
      if (line.empty()) continue;
      setProfile(cp, parseInts(line), tr);
      printf("%u\n", cp.calcXquantile(q));
    }
    return 0;
  }

  if (mode == "dropsnew") {
    KmerTranslator tr(argv[2]);
    double threshold = atof(argv[3]);
    MapTable table;
    CountProfile cp(&tr, &table);
    while (std::getline(std::cin, line)) {
      if (line.empty()) continue;
      setProfile(cp, parseInts(line), tr);
      printf("%d\n", cp.checkForSpuriousTransitionDropsWithWindowNew(threshold) ? 1 : 0);
    }
    return 0;
  }

  if (mode == "dropswindow" || mode == "dropsold") {
    KmerTranslator tr(argv[2]);
    MapTable table;
    CountProfile cp(&tr, &table);
    while (std::getline(std::cin, line)) {
      if (line.empty()) continue;
      setProfile(cp, parseInts(line), tr);
      uint32_t *mx = cp.maximize();
      int res;
      if (mode == "dropswindow")
        res = cp.checkForSpuriousTransitionDropsWithWindow(mx, atoi(argv[3]), atof(argv[4]),
                                                           atof(argv[5]), atoi(argv[6]) != 0);
      else
        res = cp.checkForSpuriousTransitionDrops(mx, atoi(argv[3]), atoi(argv[4]) != 0);
      printf("%d\n", res ? 1 : 0);
      delete[] mx;
    }
    return 0;
  }

  /* Full correction on a controlled table built from a reference sequence.
     stdin: <reference> <read> <qual|-> <multiplicity> per line. */
  if (mode == "correct" || mode == "indel" || mode == "trim") {
    KmerTranslator tr(argv[2]);
    double threshold = atof(argv[3]);
    unsigned int pseudo = atoi(argv[4]), lower = atoi(argv[5]);
    int flag = atoi(argv[6]);
    while (std::getline(std::cin, line)) {
      if (line.empty()) continue;
      std::istringstream is(line);
      std::string ref, read, qual;
      unsigned int mult;
      is >> ref >> read >> qual >> mult;
      if (qual == "-") qual = "";

      MapTable table;
      buildTableFromReference(table, tr, ref, mult);
      SequenceInfo si{"t", "", read, qual, qual.empty() ? '>' : '@'};
      CountProfile cp(&tr, &table);
      if (read.size() < tr.getSpan()) { printf("%s\t%s\t0\t0\t0\n", read.c_str(), qual.c_str()); continue; }
      cp.fill(&si);
      uint32_t *mx = cp.maximize();
      unsigned int sub = 0, ins = 0, del = 0, trimmed = 0;
      int status = ERROR_FREE;
      if (mode == "correct") {
        status = cp.doSubstitutionCorrection(mx, threshold, pseudo, lower, flag != 0, false, &sub);
        printf("%s\t%s\t%u\t%u\t%u\t%d\n", si.seq.c_str(),
               si.qual.empty() ? "-" : si.qual.c_str(), sub, ins, del, status);
      } else if (mode == "indel") {
        cp.doIndelCorrection(mx, threshold, pseudo, lower, flag != 0, false, &sub, &ins, &del);
        printf("%s\t%s\t%u\t%u\t%u\n", si.seq.c_str(),
               si.qual.empty() ? "-" : si.qual.c_str(), sub, ins, del);
      } else {
        cp.doTrimming(mx, threshold, pseudo, lower, (unsigned int) flag, false, &trimmed);
        printf("%s\t%s\t%u\n", si.seq.c_str(),
               si.qual.empty() ? "-" : si.qual.c_str(), trimmed);
      }
      delete[] mx;
    }
    return 0;
  }

  if (mode == "filename" || mode == "fileext") {
    while (std::getline(std::cin, line)) {
      /* an empty input line is a legitimate case for both */
      if (mode == "filename")
        printf("%s\n", getFilename(line).c_str());
      else
        printf("%s\n", getFileExtension(line).c_str());
    }
    return 0;
  }

  if (mode == "endswith") {
    while (std::getline(std::cin, line)) {
      size_t sp = line.find(' ');
      std::string suffix = line.substr(0, sp);
      std::string str = sp == std::string::npos ? std::string("") : line.substr(sp + 1);
      if (suffix == "@") suffix = "";
      if (str == "@") str = "";
      printf("%d\n", endsWith(suffix, str) ? 1 : 0);
    }
    return 0;
  }

  fprintf(stderr, "unknown mode %s\n", mode.c_str());
  return 2;
}
