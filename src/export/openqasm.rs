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

use error;

/// Trait for gates that can be represented in OpenQasm.
pub trait OpenQasm
{
    /// OpenQasm representation
    ///
    /// Return an OpenQasm instruction string for this gate operating on qubits
    /// `bits`. The array `bit_names` contains the names of all qubits.
    fn open_qasm(&self, bit_names: &[String], bits: &[usize]) -> String;

    /// OpenQasm representation of conditional gate.
    ///
    /// Return the OpenQasm representation of a gate that is only executed when
    /// the condition `condition` on the classical bits of the program state
    /// holds. The default implementation only works for a single gate,
    /// composite gates (like `Composite` or `Kron`) should overwrite this
    /// default. On success, returns `Ok` with the instruction string. On error,
    /// returns `Err` with an error message.
    fn conditional_open_qasm(&self, condition: &str, bit_names: &[String],
        bits: &[usize]) -> error::Result<String>
    {
        Ok(format!("if ({}) {}", condition, self.open_qasm(bit_names, bits)))
    }
}

#[cfg(test)]
mod tests
{
    use gates::H;

    use super::OpenQasm;

    #[test]
    fn test_conditional_open_qasm()
    {
        let bit_names = [String::from("qb0"), String::from("qb1")];

        let res = H::new().conditional_open_qasm("b == 0", &bit_names, &[1]);
        assert_eq!(res, Ok(String::from("if (b == 0) h qb1")));
    }
}
