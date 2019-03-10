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

/// Enumeration for errors in the use of q1tsim
#[derive(Debug, PartialEq)]
pub enum Error
{
    /// Invalid index for quantum bit
    InvalidQBit(usize),
    /// Invalid index for classical bit
    InvalidCBit(usize),
    /// Results asked for circuit that has not been run yet
    NotExecuted,
    /// Not using a full register where one is exepected
    IncompleteRegister,
    /// Trying to peek into the quantum state where it's not possible (export, mainly)
    PeekInvalid(&'static str),
    /// Classical registers not available in c-Qasm
    NoClassicalRegister,
    /// Other errors that should not occur
    InternalError(String)
}

impl ::std::fmt::Display for Error
{
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result
    {
        match *self
        {
            Error::InvalidQBit(bit) => {
                write!(f, "Invalid index {} for a quantum bit", bit)
            },
            Error::InvalidCBit(bit) => {
                write!(f, "Invalid index {} for a classical bit", bit)
            },
            Error::NotExecuted => {
                write!(f, "The circuit has not been executed yet")
            },
            Error::IncompleteRegister => {
                write!(f, "OpenQasm can only perform conditional operations based on a complete classical register")
            },
            Error::PeekInvalid(method) => {
                write!(f, "Peeking into the quantum state is not a physical operation, and is not supported in {}", method)
            }
            Error::NoClassicalRegister => {
                write!(f, "In cQasm, no classical registers can be specified. Measurements must be made to a classical bit with the same index as the qubit")
            },
            Error::InternalError(ref err) => {
                write!(f, "Internal error: {}", err)
            }
        }
    }
}

/// Type alias for a result with a q1tsim error
pub type Result<T> = ::std::result::Result<T, Error>;
