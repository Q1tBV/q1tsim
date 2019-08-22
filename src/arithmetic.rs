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
