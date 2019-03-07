// Copyright 2019 Q1t BV
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

//! The default HashMap in Rust uses a random hasher for protection against DoS
//! attacks. While this may be a good default, it is somewhat slower, especially
//! for the case of hash maps where a `u64` or `usize` is used as the key. If we
//! are not worried about malicious users, we can simply the provided key itself
//! as the hash value. Since the hash keys used in `q1tsim` are obtained from
//! random sampling, we do exactly that this module, and provide type aliases
//! for maps using the identity hasher, and with u64 or usize as value.

/// Identity hasher
///
/// Identity hasher that simply returns the value it is provided with.
pub struct IdentityHasher
{
    /// The collected hash value
    hash: u64
}

impl IdentityHasher
{
    /// Create a new identity hasher, with has value 0.
    fn new() -> Self
    {
        IdentityHasher { hash: 0 }
    }
}

impl ::std::hash::Hasher for IdentityHasher
{
    fn finish(&self) -> u64
    {
        self.hash
    }

    /// General write function
    ///
    /// This function is required by the Hasher interface, but remains
    /// unimplemented here as we only support `u64` or `usize` keys. For those
    /// types, `write_u64()` and `write_usize()` are overridden from the
    /// default.
    fn write(&mut self, _bytes: &[u8])
    {
        unimplemented!();
    }

    /// Set the hash value to `x`.
    fn write_u64(&mut self, x: u64)
    {
        self.hash = x;
    }

    /// Set the hash value to `x`.
    fn write_usize(&mut self, x: usize)
    {
        self.hash = x as u64;
    }
}

/// Structure to build an [IdentityHasher](struct.IdentityHasher.html).
pub struct BuildIdentityHasher {}

impl BuildIdentityHasher
{
    /// Create a new BuildIdentityHasher
    fn new() -> Self
    {
        BuildIdentityHasher {}
    }
}

impl ::std::hash::BuildHasher for BuildIdentityHasher
{
    type Hasher = IdentityHasher;
    fn build_hasher(&self) -> Self::Hasher
    {
        IdentityHasher::new()
    }
}

/// Type alias for HashMap using `u64` as keys, and `IdentityHasher` as hash algorithm.
pub type U64HashMap<V> = ::std::collections::HashMap<u64, V, BuildIdentityHasher>;
/// Type alias for HashMap using `usize` as keys, and `IdentityHasher` as hash algorithm.
pub type USizeHashMap<V> = ::std::collections::HashMap<usize, V, BuildIdentityHasher>;

/// Create a new and empty hash map using `u64` as keys, and `IdentityHasher` as
/// hash algorithm.
pub fn new_u64_hash_map<V>() -> U64HashMap<V>
{
    U64HashMap::with_hasher(BuildIdentityHasher::new())
}

/// Create a new and empty hash map using `usize` as keys, and `IdentityHasher` as
/// hash algorithm.
pub fn new_usize_hash_map<V>() -> USizeHashMap<V>
{
    USizeHashMap::with_hasher(BuildIdentityHasher::new())
}

#[cfg(test)]
mod tests
{
    use super::{new_u64_hash_map, new_usize_hash_map};

    #[test]
    fn test_u64()
    {
        let mut map = new_u64_hash_map();
        map.insert(1, "one");
        map.insert(2, "two");
        map.insert(0xffffffffff, "lots");

        assert_eq!(map.len(), 3);
        assert_eq!(map.keys().collect::<Vec<&u64>>(), vec![&1, &2, &0xffffffffff]);
        assert_eq!(map[&1], "one");
        assert_eq!(map[&2], "two");
        assert_eq!(map[&0xffffffffff], "lots");
    }

    #[test]
    fn test_usize()
    {
        let mut map = new_usize_hash_map();
        map.insert(1, "one");
        map.insert(2, "two");
        map.insert(0xffffffff, "many");

        assert_eq!(map.len(), 3);
        assert_eq!(map.keys().collect::<Vec<&usize>>(), vec![&1, &2, &0xffffffff]);
        assert_eq!(map[&1], "one");
        assert_eq!(map[&2], "two");
        assert_eq!(map[&0xffffffff], "many");
    }
}
