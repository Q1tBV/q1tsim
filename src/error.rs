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
    /// Unable to create condition operation in c-Qasm
    InvalidConditionalOp(String),
    /// Trying to export in a format for which no export function was written
    NotImplemented(&'static str, String),
    /// Trying to close a loop where none is open in LaTeX export
    CantCloseLoop,
    /// Trying to reserve range in LaTeX export, but previous reservation is open
    RangeAlreadyOpen
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
            ExportError::InvalidConditionalOp(ref qasm) => {
                write!(f, "Unable to find gate name or argument in \"{}\"", qasm)
            },
            ExportError::NotImplemented(method, ref desc) => {
                write!(f, "Export to {} was not implemented for \"{}\"", method, desc)
            },
            ExportError::CantCloseLoop => {
                write!(f, "Unable to close loop, because no loop is currently open")
            },
            ExportError::RangeAlreadyOpen => {
                write!(f, "Trying to reserve range of bits, but a previous reservation is still open")
            }
        }
    }
}

/// Type alias for a result with an export error
pub type ExportResult<T> = ::std::result::Result<T, ExportError>;

/// Structure for errors encountered while parsing a composite gate description
#[derive(Debug, PartialEq)]
pub enum ParseError
{
    /// Gate name not recognised
    UnknownGate(String),
    /// No gate name found
    NoGateName(String),
    /// Wrong number of arguments to gate
    InvalidNrArguments(usize, usize, String),
    /// Invalid number of qubits to operate on
    InvalidNrBits(usize, usize, String),
    /// Unable to parse argument to gate
    InvalidArgument(String),
    /// Unable to find bit numbers on which the gate operates
    NoBits(String),
    /// Unable to parse bit number
    InvalidBit(String),
    /// Text occurs after a gate description
    TrailingText(String),
    /// Unclosed parentheses in argument expression
    UnclosedParentheses(String),
}

impl ::std::fmt::Display for ParseError
{
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result
    {
        match *self
        {
            ParseError::UnknownGate(ref name) => {
                write!(f, "Unknown gate \"{}\"", name)
            },
            ParseError::NoGateName(ref text) => {
                write!(f, "Failed to find gate name in \"{}\"", text)
            },
            ParseError::InvalidNrArguments(actual, expected, ref name) => {
                write!(f, "Expected {} arguments to \"{}\" gate, got {}", expected, name, actual)
            },
            ParseError::InvalidNrBits(actual, expected, ref name) => {
                write!(f, "Expected {} bits for \"{}\" gate, got {}", expected, name, actual)
            },
            ParseError::InvalidArgument(ref text) => {
                write!(f, "Failed to parse argument \"{}\"", text)
            },
            ParseError::NoBits(ref name) => {
                write!(f, "Unable to find the bits gate {} operates on", name)
            },
            ParseError::InvalidBit(ref text) => {
                write!(f, "Failed to parse bit number in \"{}\"", text)
            },
            ParseError::TrailingText(ref text) => {
                write!(f, "Trailing text after gate description: \"{}\"", text)
            },
            ParseError::UnclosedParentheses(ref text) => {
                write!(f, "Unclosed parentheses in expression: \"{}\"", text)
            }
        }
    }
}

/// Type alias for a result with a parse error
pub type ParseResult<T> = ::std::result::Result<T, ParseError>;

/// Enumeration for errors in the use of q1tsim
#[derive(Debug, PartialEq)]
pub enum Error
{
    /// Number of bits passed does not match gate
    InvalidNrBits(usize, usize, String),
    /// Invalid index for quantum bit
    InvalidQBit(usize),
    /// Invalid index for classical bit
    InvalidCBit(usize),
    /// Using classically controlled operation with wrong number of control bits
    InvalidNrControlBits(usize, usize, String),
    /// Trying to store measurements into arry that cannot hold them
    NotEnoughSpace(usize, usize),
    /// Results asked for circuit that has not been run yet
    NotExecuted,
    /// Acting with a non-stabilizer gate on a stabilizer circuit
    NotAStabilizer(String),
    /// Other errors that should not occur
    InternalError(String),
    /// Error reating to the export of a circuit
    ExportError(ExportError),
    /// Error in parsing a composite gate description
    ParseError(ParseError)
}

impl From<ExportError> for Error
{
    fn from(err: ExportError) -> Self
    {
        Error::ExportError(err)
    }
}

impl From<ParseError> for Error
{
    fn from(err: ParseError) -> Self
    {
        Error::ParseError(err)
    }
}

impl ::std::fmt::Display for Error
{
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result
    {
        match *self
        {
            Error::InvalidNrBits(actual, expected, ref desc) => {
                write!(f, "Expected {} bits for \"{}\", got {}", expected, desc, actual)
            },
            Error::InvalidQBit(bit) => {
                write!(f, "Invalid index {} for a quantum bit", bit)
            },
            Error::InvalidCBit(bit) => {
                write!(f, "Invalid index {} for a classical bit", bit)
            },
            Error::InvalidNrControlBits(actual, expected, ref desc) => {
                write!(f, "The number of runs is {}, but received {} control bits for controlled {} operation",
                    expected, actual, desc)
            },
            Error::NotEnoughSpace(actual, expected) => {
                write!(f, "Not enough space to store {} measurement results in array of length {}",
                    expected, actual)
            },
            Error::NotExecuted => {
                write!(f, "The circuit has not been executed yet")
            },
            Error::NotAStabilizer(ref desc) => {
                write!(f, "{} is no a stabilizer gate", desc)
            },
            Error::InternalError(ref err) => {
                write!(f, "Internal error: {}", err)
            },
            Error::ExportError(ref err) => {
                write!(f, "{}", err)
            },
            Error::ParseError(ref err) => {
                write!(f, "{}", err)
            }
        }
    }
}

/// Type alias for a result with a q1tsim error
pub type Result<T> = ::std::result::Result<T, Error>;
