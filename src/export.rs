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

use crate::gates::Gate;

mod cqasm;
mod latex;
mod openqasm;

pub use self::cqasm::CQasm;
pub use self::latex::{Latex, LatexExportState};
pub use self::openqasm::OpenQasm;

/// Trait combining the traits necessary for a gate in a quantum circuit
pub trait CircuitGate: Gate + OpenQasm + CQasm + Latex
{
    fn as_gate(&self) -> &dyn Gate;
    fn clone_box(&self) -> Box<dyn CircuitGate>;
}

impl<G: 'static + Clone + Gate + OpenQasm + CQasm + Latex> CircuitGate for G
{
    fn as_gate(&self) -> &dyn Gate
    {
        self
    }

    fn clone_box(&self) -> Box<dyn CircuitGate>
    {
        Box::new(self.clone())
    }
}

impl Clone for Box<dyn CircuitGate>
{
    fn clone(&self) -> Self
    {
        self.clone_box()
    }
}
