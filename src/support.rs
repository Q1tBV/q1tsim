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
