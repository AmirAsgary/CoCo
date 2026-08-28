// Written by Annika Jochheim <annika.jochheim@mpinat.mpg.de>

#include "KmerTranslator.h"
#include "Info.h"
#include "util.h"

KmerTranslator::KmerTranslator(std::string spacedKmerPattern) {

  span = spacedKmerPattern.length();
  if (span > 64) {
    Info(Info::INFO) << "Error: Only kmerSpan <= 64 is supported\n";
    EXIT(EXIT_FAILURE);
  }
  weight = 0;
  for (unsigned int idx = 0; idx < span; idx++) {
    if (spacedKmerPattern[idx] == '1') {
      weight++;
    } else if (spacedKmerPattern[idx] == '0')
      continue;
    else {
      Info(Info::INFO) << "Error: Found invalid character in spacedKmerPattern\n";
      EXIT(EXIT_FAILURE);
    }
    if (spacedKmerPattern[idx] != spacedKmerPattern[span - idx - 1]) {
      Info(Info::INFO) << "Error: spacedKmerPattern must be symmetric\n";
      EXIT(EXIT_FAILURE);
    }
  }

  if (spacedKmerPattern[0] != '1' || spacedKmerPattern[span - 1] != '1') {
    Info(Info::INFO) << "Error: First and last position in spacedKmerPattern must 1\n";
    EXIT(EXIT_FAILURE);
  }

  if (weight < 12 || weight > 32) {
    Info(Info::INFO) << "Error: Only 12 <= weight <= 32 is supported\n";
    EXIT(EXIT_FAILURE);
  }

  _mask_array = new unsigned char[weight];
  _inverse_mask_array = new unsigned char[span];
  unsigned int jdx = 0;
  for (unsigned int idx = 0; idx < span; idx++) {
    if (spacedKmerPattern[idx] == '1') {
      _mask_array[jdx] = idx;
      _inverse_mask_array[idx] = jdx;
      jdx++;
    } else {
      _inverse_mask_array[idx] = UCHAR_MAX;
    }
  }

  _spaced_mask = std::stoll(spacedKmerPattern, nullptr, 2);

}

KmerTranslator::~KmerTranslator() {
  delete[] _mask_array;
  delete[] _inverse_mask_array;
}

unsigned short KmerTranslator::getSpan() const {
  return span;
}

unsigned short KmerTranslator::getWeight() const {
  return weight;
}

unsigned short KmerTranslator::getLongestBlock() const {

  unsigned short blockCounter = 0, maxBlockCounter = 0;

  for (size_t idx = 0; idx < this->span; idx++) {
    if(_inverse_mask_array[idx]!= UCHAR_MAX)
      blockCounter++;
    else
      blockCounter = 0;

    if(blockCounter > maxBlockCounter)
      maxBlockCounter = blockCounter;
  }
  return maxBlockCounter;
}

packedKmerType KmerTranslator::kmer2packedKmer(const spacedKmerType kmer) const {

  packedKmerType packedKmer = 0;

  /* PERF: the original walked all `span` positions and tested `i == _mask_array[j]`
   * on each. But _mask_array already lists exactly the `weight` informative offsets,
   * in ascending order (see the constructor), so we can walk it directly. That drops
   * span-weight iterations (33->27 for the unik mask, 41->32 for CoCo's default) and,
   * more importantly, removes a data-dependent branch that mispredicts in an
   * irregular pattern set by the mask. The bit extraction is also simplified from
   * (kmer & (3 << shift)) >> shift to (kmer >> shift) & 3, which is equivalent.
   * The emitted bit order is unchanged, so packed k-mers are identical. */
  for (size_t j = 0; j < weight; j++) {
    const unsigned int shift = 2u * (unsigned int)(span - 1 - _mask_array[j]);
    packedKmer = (packedKmer << 2) | (packedKmerType) ((kmer >> shift) & (spacedKmerType) 3);
  }

  return packedKmer;
}

packedKmerType KmerTranslator::kmer2minPackedKmer(const spacedKmerType kmer) const {
  packedKmerType packedKmer = kmer2packedKmer(kmer);
  return (minIndex(packedKmer, weight));
}

packedKmerType KmerTranslator::kmer2minPackedKmer(const packedKmerType kmer) const {
  return (minIndex(kmer, weight));
}

void KmerTranslator::getBestSplit(unsigned int &logIndexSize, unsigned int &logOffsetSize ) const {

  if (2*weight > 30) {
    logIndexSize = 30;
    logOffsetSize = 2 * weight - logIndexSize;
  } else {
    logIndexSize = 2 * weight-2;
    logOffsetSize = 2;
  }

}