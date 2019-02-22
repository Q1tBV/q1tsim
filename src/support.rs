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

/// Combine indexes into ranges of indexes
///
/// Combine the index in `nrs` into ranges. Each range is a pair of indices
/// `(first, last)`, with each number in the range `[first..last]` (inclusive)
/// present in `nrs`.
pub fn get_ranges(nrs: &[usize]) -> Vec<(usize, usize)>
{
    let mut snrs = nrs.to_vec();
    snrs.sort();

    let mut first = snrs[0];
    let mut last = first;
    let mut ranges = vec![];
    for &nr in snrs[1..].iter()
    {
        if nr == last+1
        {
            last += 1;
        }
        else
        {
            ranges.push((first, last));
            first = nr;
            last = first;
        }
    }
    ranges.push((first, last));

    ranges
}

#[cfg(test)]
mod tests
{
    use super::get_ranges;

    #[test]
    fn test_get_ranges()
    {
        let ranges = get_ranges(&[0, 1, 2, 3, 4, 5, 6]);
        assert_eq!(ranges, vec![(0, 6)]);
        let ranges = get_ranges(&[0, 1, 3, 4, 5, 6]);
        assert_eq!(ranges, vec![(0, 1), (3, 6)]);
        let ranges = get_ranges(&[0, 2, 3, 4, 5, 6]);
        assert_eq!(ranges, vec![(0, 0), (2, 6)]);
        let ranges = get_ranges(&[0, 1, 2, 3, 4, 6]);
        assert_eq!(ranges, vec![(0, 4), (6, 6)]);
        let ranges = get_ranges(&[2, 4, 3, 1, 0, 6]);
        assert_eq!(ranges, vec![(0, 4), (6, 6)]);
        let ranges = get_ranges(&[5, 6, 9, 2, 0, 11]);
        assert_eq!(ranges, vec![(0, 0), (2, 2), (5, 6), (9, 9), (11, 11)]);
    }
}
