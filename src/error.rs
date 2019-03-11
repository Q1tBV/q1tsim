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

/// Enumeration for errors relating to the export of circuits
#[derive(Debug, PartialEq)]
pub enum ExportError
{
    /// Trying to export a peek into the quantum state
    ExportPeekInvalid(&'static str),
    /// Classical registers not available in c-Qasm
    NoClassicalRegister,
    /// Not using a full condition register in OpenQasm export
    IncompleteConditionRegister,
}

impl ::std::fmt::Display for ExportError
{
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result
    {
        match *self
        {
            ExportError::ExportPeekInvalid(method) => {
                write!(f, "Peeking into the quantum state is not a physical operation, and is not supported in {}", method)
            },
            ExportError::NoClassicalRegister => {
                write!(f, "In cQasm, no classical registers can be specified. Measurements must be made to a classical bit with the same index as the qubit")
            },
            ExportError::IncompleteConditionRegister => {
                write!(f, "OpenQasm can only perform conditional operations based on a complete classical register")
            },
        }
    }
}

/// Type alias for a result with an export error
pub type ExportResult<T> = ::std::result::Result<T, ExportError>;

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
    /// Other errors that should not occur
    InternalError(String),
    /// Error reating to the export of a circuit
    ExportError(ExportError)
}

impl From<ExportError> for Error
{
    fn from(err: ExportError) -> Self
    {
        Error::ExportError(err)
    }
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
            Error::InternalError(ref err) => {
                write!(f, "Internal error: {}", err)
            },
            Error::ExportError(ref err) => {
                write!(f, "{}", err)
            }
        }
    }
}

/// Type alias for a result with a q1tsim error
pub type Result<T> = ::std::result::Result<T, Error>;
