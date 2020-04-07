// Copyright (C) 2020 ETH Zurich
// Copyright (C) 2020 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Provides a set with O(log n) insertion, removal and random access, i.e. access to the i-th lowest key.
// Useful for a random selection in an ordered list with variable size.
// Implemented as an augmented red-black tree.
// Code to fix RB color violation from https://www.geeksforgeeks.org/red-black-tree-set-1-introduction-2/

#ifndef DCA_UTIL_TYPE_RANDOM_ACCESS_SET_HPP
#define DCA_UTIL_TYPE_RANDOM_ACCESS_SET_HPP

#include "dca/util/containers/random_access_map.hpp"

namespace dca {
namespace util {

// Precondition: elements of type Key have full order.
template <class Key, std::size_t chunk_size = 256>
class RandomAccessSet {
public:
  RandomAccessSet() = default;
  RandomAccessSet(const RandomAccessSet& rhs) = default;
  RandomAccessSet(RandomAccessSet&& rhs) = default;
  RandomAccessSet& operator=(const RandomAccessSet& rhs) = default;
  RandomAccessSet& operator=(RandomAccessSet&& rhs) = default;

  RandomAccessSet(const std::initializer_list<Key>& list);

  // Insert new key.
  // Precondition: key is not already in the map.
  void insert(const Key& key);

  // Remove the node relative to key.
  // Precondition: key is in the map.
  void erase(const Key& key);

  // Returns a reference to the value relative to the i-th key.
  // Precondition: 0 <= index < size()
  const Key& operator[](const std::size_t index) const;

  // Number of keys stored in the map.
  std::size_t size() const {
    return map_.size();
  }

  // Returns an array of ordered keys and value pairs.
  std::vector<Key> linearize() const;

  bool checkConsistency() const {
    return map_.checkConsistency();
  }

private:
  struct Null {
    Null() = default;
  };

  // Members
  dca::util::RandomAccessMap<Key, Null, chunk_size> map_;
};

template <class Key, std::size_t chunk_size>
RandomAccessSet<Key, chunk_size>::RandomAccessSet(const std::initializer_list<Key>& list) {
  for (const auto k : list)
    map_.insert(k, {});
}

template <class Key, std::size_t chunk_size>
void RandomAccessSet<Key, chunk_size>::insert(const Key& key) {
  map_.insert(key, {});
}

template <class Key, std::size_t chunk_size>
void RandomAccessSet<Key, chunk_size>::erase(const Key& key) {
  map_.erase(key);
}

template <class Key, std::size_t chunk_size>
const Key& RandomAccessSet<Key, chunk_size>::operator[](const std::size_t index) const {
  const auto& [key, _] = map_[index];
  return key;
}

template <class Key, std::size_t chunk_size>
std::vector<Key> RandomAccessSet<Key, chunk_size>::linearize() const {
  auto linear_map = map_.linearize();
  std::vector<Key> result;
  result.reserve(linear_map.size());
  for (const auto& [key, _] : linear_map)
    result.push_back(key);

  return result;
}

}  // namespace util
}  // namespace dca

#endif  // #define DCA_UTIL_TYPE_RANDOM_ACCESS_SET_HPP
