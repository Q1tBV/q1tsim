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

/// Trait for squaring a gate
///
/// Some quantum algorithms use the square or higher powers of a gate. If
/// the square of a gate can be represted more simply than applying the
/// gate twice, gates can override this function. The default returns an
/// `OpNotImplemented` error.
pub trait Square: crate::gates::Gate
{
    /// The result type of the square of this gate
    type SqType;

    /// Square of this gate
    fn square(&self) -> crate::error::Result<Self::SqType>
    {
        Err(crate::error::Error::OpNotImplemented(String::from("square"),
            String::from(self.description())))
    }
}
