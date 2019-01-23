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


use gates;

use gates::Gate;

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
        bits: &[usize]) -> Result<String, String>
    {
        Ok(format!("if ({}) {}", condition, self.open_qasm(bit_names, bits)))
    }
}

/// Trait for gates that can be represented in c-Qasm.
pub trait CQasm
{
    /// cQasm representation
    ///
    /// Return an cQasm instruction string for this gate operating on qubits
    /// `bits`. The array `bit_names` contains the names of all qubits.
    fn c_qasm(&self, bit_names: &[String], bits: &[usize]) -> String;

    /// cQasm representation of conditional gate.
    ///
    /// Return the cQasm representation of a gate that is only executed when
    /// the condition `condition` on the classical bits of the program state
    /// holds. The default implementation only works for a single gate,
    /// composite gates (like `Composite` or `Kron`) should overwrite this
    /// default. On success, returns `Ok` with the instruction string. On error,
    /// returns `Err` with an error message.
    fn conditional_c_qasm(&self, condition: &str, bit_names: &[String],
        bits: &[usize]) -> Result<String, String>
    {
        let unc_qasm = self.c_qasm(bit_names, bits);
        let parts: Vec<&str> = unc_qasm.splitn(2, " ").collect();
        if parts.len() != 2
        {
            Err(format!("Unable to find instruction name or argument in \"{}\"", unc_qasm))
        }
        else
        {
            Ok(format!("c-{} {} {}", parts[0], condition, parts[1]))
        }
    }
}


/// Trait combining the traits necessary for a gate in a quantum circuit
pub trait CircuitGate: gates::Gate + OpenQasm + CQasm {}

impl<G: Gate + OpenQasm + CQasm> CircuitGate for G {}
