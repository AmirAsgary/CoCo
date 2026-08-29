//! Monotonic-deque sliding maximum.
//!
//! `calcNeighborhoodTolerance` keeps a ~`span`-wide window over the maximized
//! profile and calls `std::max_element` on it once per position -- for CoCo's
//! default pattern that is a 42-element linear scan at each of ~150 positions, for
//! every read, and the tolerance array is rebuilt on every correction round.
//!
//! This gives the same answer in O(1) amortised per operation. It supports an
//! arbitrary interleaving of `push` and `pop_front`, which matters because the
//! original does not maintain a clean fixed-width window: for a stretch of
//! positions it pushes a value it has already pushed, leaving a duplicate in the
//! window. Sequence numbers rather than positions identify the logical front, so
//! that behaviour is reproduced rather than tidied up.

use std::collections::VecDeque;

pub struct SlidingMax {
    /// Values in non-increasing order, tagged with insertion sequence.
    dq: VecDeque<(u32, u64)>,
    next_seq: u64,
    head_seq: u64,
}

impl SlidingMax {
    pub fn new() -> Self {
        SlidingMax { dq: VecDeque::new(), next_seq: 0, head_seq: 0 }
    }

    pub fn with_capacity(n: usize) -> Self {
        SlidingMax { dq: VecDeque::with_capacity(n), next_seq: 0, head_seq: 0 }
    }

    pub fn clear(&mut self) {
        self.dq.clear();
        self.next_seq = 0;
        self.head_seq = 0;
    }

    /// Append a value to the logical window.
    #[inline]
    pub fn push(&mut self, v: u32) {
        while let Some(&(back, _)) = self.dq.back() {
            if back <= v {
                self.dq.pop_back();
            } else {
                break;
            }
        }
        self.dq.push_back((v, self.next_seq));
        self.next_seq += 1;
    }

    /// Drop the oldest value from the logical window.
    #[inline]
    pub fn pop_front(&mut self) {
        if let Some(&(_, seq)) = self.dq.front() {
            if seq == self.head_seq {
                self.dq.pop_front();
            }
        }
        self.head_seq += 1;
    }

    /// Largest value currently in the window.
    #[inline]
    pub fn max(&self) -> u32 {
        self.dq.front().map(|&(v, _)| v).unwrap_or(0)
    }

    #[inline]
    pub fn len(&self) -> usize {
        (self.next_seq - self.head_seq) as usize
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

impl Default for SlidingMax {
    fn default() -> Self {
        Self::new()
    }
}

/// Same interface, tracking the minimum instead. Used by the filter windows.
pub struct SlidingMin {
    dq: VecDeque<(u32, u64)>,
    next_seq: u64,
    head_seq: u64,
}

impl SlidingMin {
    pub fn new() -> Self {
        SlidingMin { dq: VecDeque::new(), next_seq: 0, head_seq: 0 }
    }
    pub fn clear(&mut self) {
        self.dq.clear();
        self.next_seq = 0;
        self.head_seq = 0;
    }
    #[inline]
    pub fn push(&mut self, v: u32) {
        while let Some(&(back, _)) = self.dq.back() {
            if back >= v {
                self.dq.pop_back();
            } else {
                break;
            }
        }
        self.dq.push_back((v, self.next_seq));
        self.next_seq += 1;
    }
    #[inline]
    pub fn pop_front(&mut self) {
        if let Some(&(_, seq)) = self.dq.front() {
            if seq == self.head_seq {
                self.dq.pop_front();
            }
        }
        self.head_seq += 1;
    }
    #[inline]
    pub fn min(&self) -> u32 {
        self.dq.front().map(|&(v, _)| v).unwrap_or(0)
    }
    #[inline]
    pub fn len(&self) -> usize {
        (self.next_seq - self.head_seq) as usize
    }
}

impl Default for SlidingMin {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Reference: a plain Vec with `max_element`, i.e. exactly what the C++ does.
    struct Naive(Vec<u32>);
    impl Naive {
        fn push(&mut self, v: u32) {
            self.0.push(v);
        }
        fn pop_front(&mut self) {
            self.0.remove(0);
        }
        fn max(&self) -> u32 {
            self.0.iter().copied().max().unwrap_or(0)
        }
        fn min(&self) -> u32 {
            self.0.iter().copied().min().unwrap_or(0)
        }
    }

    #[test]
    fn matches_naive_under_random_operations() {
        let mut state: u64 = 0x1234_5678_9ABC_DEF0;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        for _ in 0..200 {
            let mut sm = SlidingMax::new();
            let mut naive = Naive(Vec::new());
            for _ in 0..500 {
                let op = next() % 3;
                if op == 0 && !naive.0.is_empty() {
                    sm.pop_front();
                    naive.pop_front();
                } else {
                    let v = (next() % 1000) as u32;
                    sm.push(v);
                    naive.push(v);
                }
                assert_eq!(sm.len(), naive.0.len());
                assert_eq!(sm.max(), naive.max());
            }
        }
    }

    #[test]
    fn min_variant_matches_naive() {
        let mut state: u64 = 0xFEED_FACE_CAFE_BEEF;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        let mut sm = SlidingMin::new();
        let mut naive = Naive(Vec::new());
        for _ in 0..5000 {
            if next() % 3 == 0 && !naive.0.is_empty() {
                sm.pop_front();
                naive.pop_front();
            } else {
                let v = (next() % 100) as u32;
                sm.push(v);
                naive.push(v);
            }
            assert_eq!(sm.min(), naive.min());
        }
    }

    #[test]
    fn reproduces_the_duplicate_push_pattern() {
        // The tolerance window pushes maxProfile[span/2] repeatedly at the start,
        // so the same value can sit in the window twice. Popping once must not
        // remove both copies.
        let mut sm = SlidingMax::new();
        sm.push(5);
        sm.push(5);
        sm.push(3);
        assert_eq!(sm.max(), 5);
        sm.pop_front();
        assert_eq!(sm.max(), 5);
        sm.pop_front();
        assert_eq!(sm.max(), 3);
    }
}
