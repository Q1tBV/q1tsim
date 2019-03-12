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
use gates;

/// Trait for gates that can be represented in c-Qasm.
pub trait CQasm: gates::Gate
{
    /// cQasm representation
    ///
    /// Return an cQasm instruction string for this gate operating on qubits
    /// `bits`. The array `bit_names` contains the names of all qubits. The
    /// default implementation returns a NotImplemented error.
    fn c_qasm(&self, _bit_names: &[String], _bits: &[usize])
        -> error::ExportResult<String>
    {
        Err(error::ExportError::NotImplemented("c-Qasm",
            String::from(self.description())))
    }

    /// cQasm representation of conditional gate.
    ///
    /// Return the cQasm representation of a gate that is only executed when
    /// the condition `condition` on the classical bits of the program state
    /// holds. The default implementation only works for a single gate,
    /// composite gates (like `Composite` or `Kron`) should overwrite this
    /// default. On success, returns `Ok` with the instruction string. On error,
    /// returns `Err` with an error message.
    fn conditional_c_qasm(&self, condition: &str, bit_names: &[String],
        bits: &[usize]) -> error::ExportResult<String>
    {
        let unc_qasm = self.c_qasm(bit_names, bits)?;
        let parts: Vec<&str> = unc_qasm.splitn(2, " ").collect();
        if parts.len() != 2
        {
            // This shouldn't happen, really.
            Err(error::ExportError::InvalidConditionalOp(unc_qasm.clone()))
        }
        else
        {
            Ok(format!("c-{} {}, {}", parts[0], condition, parts[1]))
        }
    }
}

#[cfg(test)]
mod tests
{
    use gates::H;

    use super::CQasm;

    #[test]
    fn test_conditional_c_qasm()
    {
        let bit_names = [String::from("qb0"), String::from("qb1")];

        let res = H::new().conditional_c_qasm("b[0]", &bit_names, &[1]);
        assert_eq!(res, Ok(String::from("c-h b[0], qb1")));
    }
}
