use crate::circuit::{Basis, Circuit};
use crate::gates::*;
use ::std::os::raw::c_char;

const RESULT_ERROR: u32 = 0;
const RESULT_EMPTY: u32 = 1;
const RESULT_STRING: u32 = 2;
const RESULT_HISTOGRAM: u32 = 3;
const RESULT_CSTATE: u32 = 5;

fn to_cstring(msg: &str) -> *const c_char
{
    let cstring = ::std::ffi::CString::new(msg).unwrap();
    let ptr = cstring.as_ptr();
    // Prevent Rust from cleaning up this string
    ::std::mem::forget(cstring);
    ptr
}

fn cstring_free(ptr: *mut c_char)
{
    unsafe { ::std::ffi::CString::from_raw(ptr); }
}

#[repr(C)]
pub struct CHistElem
{
    key: *const c_char,
    count: usize
}

impl CHistElem
{
    fn new(key: &str, count: usize) -> Self
    {
        CHistElem { key: to_cstring(key), count: count }
    }
}

#[derive(Clone, Copy)]
#[repr(C)]
pub struct CParameter
{
    value: f64,
    value_ptr: *const f64
}

impl From<CParameter> for crate::gates::Parameter
{
    fn from(p: CParameter) -> crate::gates::Parameter
    {
        if p.value_ptr.is_null()
        {
            crate::gates::Parameter::Direct(p.value)
        }
        else
        {
            crate::gates::Parameter::FFIRef(p.value_ptr)
        }
    }
}

#[repr(C)]
pub struct CResult
{
    data: *const ::std::ffi::c_void,
    length: usize,
    size: usize,
    restype: u32
}

impl CResult
{
    fn new() -> Self
    {
        CResult
        {
            data: ::std::ptr::null(),
            length: 0,
            size: 0,
            restype: RESULT_EMPTY
        }
    }

    fn error(msg: &str) -> Self
    {
        CResult
        {
            data: to_cstring(msg) as *const ::std::ffi::c_void,
            length: 0,
            size: 0,
            restype: RESULT_ERROR
        }
    }

    fn string(msg: &str) -> Self
    {
        CResult
        {
            data: to_cstring(msg) as *const ::std::ffi::c_void,
            length: 0,
            size: 0,
            restype: RESULT_STRING
        }
    }

    fn histogram(hist: &::std::collections::HashMap<String, usize>) -> Self
    {
        let elems: Vec<CHistElem> = hist.iter().map(|(k, &v)| CHistElem::new(k, v)).collect();
        let ptr = elems.as_ptr();
        let length = elems.len();
        let size = elems.capacity();
        ::std::mem::forget(elems);
        CResult
        {
            data: ptr as *const ::std::ffi::c_void,
            length: length,
            size: size,
            restype: RESULT_HISTOGRAM
        }
    }

    fn c_state(state: &[u64]) -> Self
    {
        let elems = state.to_vec();
        let ptr = elems.as_ptr();
        let length = elems.len();
        let size = elems.capacity();
        ::std::mem::forget(elems);
        CResult
        {
            data: ptr as *const ::std::ffi::c_void,
            length: length,
            size: size,
            restype: RESULT_CSTATE
        }
    }

    fn free(self)
    {
        match self.restype
        {
            RESULT_ERROR|RESULT_STRING => {
                cstring_free(self.data as *mut c_char);
            },
            RESULT_HISTOGRAM => {
                let ptr = self.data as *mut CHistElem;
                let mut elems = unsafe { Vec::from_raw_parts(ptr, self.length, self.size) };
                for elem in elems.iter_mut()
                {
                    let ptr = std::mem::replace(&mut elem.key, ::std::ptr::null());
                    cstring_free(ptr as *mut c_char);
                }
            },
            RESULT_CSTATE => {
                let ptr = self.data as *mut u64;
                let _elems = unsafe { Vec::from_raw_parts(ptr, self.length, self.size) };
                /* and drop _elems */
            },
            RESULT_EMPTY|_ => { /* ignore */ }
        }
    }
}

#[no_mangle]
pub extern "C" fn result_free(result: CResult)
{
    result.free()
}

#[no_mangle]
pub extern "C" fn circuit_new(nr_qbits: usize, nr_cbits: usize) -> *mut Circuit
{
    let ptr = Box::new(Circuit::new(nr_qbits, nr_cbits));
    Box::into_raw(ptr)
}

#[no_mangle]
pub extern "C" fn circuit_free(ptr: *mut Circuit)
{
    if !ptr.is_null()
    {
        // Convert pointer back to box, and let it fall out of scope
        unsafe { Box::from_raw(ptr); }
    }
}

#[no_mangle]
pub extern "C" fn circuit_nr_qbits(ptr: *const Circuit) -> usize
{
    assert!(!ptr.is_null());
    let circuit = unsafe { &*ptr };
    circuit.nr_qbits()
}

#[no_mangle]
pub extern "C" fn circuit_nr_cbits(ptr: *const Circuit) -> usize
{
    assert!(!ptr.is_null());
    let circuit = unsafe { &*ptr };
    circuit.nr_cbits()
}

#[no_mangle]
pub extern "C" fn circuit_cstate(ptr: *const Circuit) -> CResult
{
    assert!(!ptr.is_null());
    let circuit = unsafe { &*ptr };
    match circuit.cstate()
    {
        Some(state) => CResult::c_state(state.as_slice().unwrap()),
        None => CResult::error("Circuit has not been run yet")
    }
}

macro_rules! add_parametrized_gate
{
    (1, $circuit:ident, $gate_type:ty, $qbits:expr, $params:ident) => {
        if $params.len() != 1
        {
            Err(crate::error::Error::from(
                crate::error::ParseError::InvalidNrArguments($params.len(), 1,
                    String::from(stringify!($gate_type))
                )
            ))
        }
        else
        {
            $circuit.add_gate(<$gate_type>::new($params[0]), $qbits)
        }
    };
    (1, $circuit:ident, $control:expr, $target:expr, $gate_type:ty, $qbits:expr, $params:ident) => {
        if $params.len() != 1
        {
            Err(crate::error::Error::from(
                crate::error::ParseError::InvalidNrArguments($params.len(), 1,
                    String::from(stringify!($gate_type))
                )
            ))
        }
        else
        {
            $circuit.add_conditional_gate($control, $target,
                <$gate_type>::new($params[0]), $qbits)
        }
    };
    (2, $circuit:ident, $gate_type:ty, $qbits:expr, $params:ident) => {
        if $params.len() != 2
        {
            Err(crate::error::Error::from(
                crate::error::ParseError::InvalidNrArguments($params.len(), 2,
                    String::from(stringify!($gate_type))
                )
            ))
        }
        else
        {
            $circuit.add_gate(<$gate_type>::new($params[0], $params[1]), $qbits)
        }
    };
    (2, $circuit:ident, $control:expr, $target:expr, $gate_type:ty, $qbits:expr, $params:ident) => {
        if $params.len() != 2
        {
            Err(crate::error::Error::from(
                crate::error::ParseError::InvalidNrArguments($params.len(), 1,
                    String::from(stringify!($gate_type))
                )
            ))
        }
        else
        {
            $circuit.add_conditional_gate($control, $target,
                <$gate_type>::new($params[0], $params[1]), $qbits)
        }
    };
    (3, $circuit:ident, $gate_type:ty, $qbits:expr, $params:ident) => {
        if $params.len() != 3
        {
            Err(crate::error::Error::from(
                crate::error::ParseError::InvalidNrArguments($params.len(), 3,
                    String::from(stringify!($gate_type))
                )
            ))
        }
        else
        {
            $circuit.add_gate(<$gate_type>::new($params[0], $params[1], $params[2]), $qbits)
        }
    };
    (3, $circuit:ident, $control:expr, $target:expr, $gate_type:ty, $qbits:expr, $params:ident) => {
        if $params.len() != 3
        {
            Err(crate::error::Error::from(
                crate::error::ParseError::InvalidNrArguments($params.len(), 1,
                    String::from(stringify!($gate_type))
                )
            ))
        }
        else
        {
            $circuit.add_conditional_gate($control, $target,
                <$gate_type>::new($params[0], $params[1], $params[2]), $qbits)
        }
    };
}

#[no_mangle]
pub extern "C" fn circuit_add_gate(ptr: *mut Circuit, gate: *const c_char,
    qbits_ptr: *const usize, nr_qbits: usize,
    param_ptr: *const CParameter, nr_params: usize) -> CResult
{
    if ptr.is_null()
    {
        CResult::error("Pointer to circuit is NULL")
    }
    else if qbits_ptr.is_null()
    {
        CResult::error("Pointer to bit indices is NULL")
    }
    else
    {
        if let Ok(gate_name) = unsafe { ::std::ffi::CStr::from_ptr(gate).to_str() }
        {
            let circuit = unsafe { &mut *ptr };
            let qbits = unsafe { ::std::slice::from_raw_parts(qbits_ptr, nr_qbits) };
            let params = if param_ptr.is_null()
                {
                    &[]
                }
                else
                {
                    unsafe { ::std::slice::from_raw_parts(param_ptr, nr_params) }
                };

            let res = match gate_name.to_lowercase().as_str()
                {
                    "ch"   => { circuit.add_gate(CH::new(), qbits) },
                    "crx"   => { add_parametrized_gate!(1, circuit, CRX, qbits, params) },
                    "cry"   => { add_parametrized_gate!(1, circuit, CRY, qbits, params) },
                    "crz"   => { add_parametrized_gate!(1, circuit, CRZ, qbits, params) },
                    "cx"   => { circuit.add_gate(CX::new(), qbits) },
                    "cy"   => { circuit.add_gate(CY::new(), qbits) },
                    "cz"   => { circuit.add_gate(CZ::new(), qbits) },
                    "h"    => { circuit.add_gate(H::new(), qbits) },
                    "i"    => { circuit.add_gate(I::new(), qbits) },
                    "rx"   => { add_parametrized_gate!(1, circuit, RX, qbits, params) },
                    "ry"   => { add_parametrized_gate!(1, circuit, RY, qbits, params) },
                    "rz"   => { add_parametrized_gate!(1, circuit, RZ, qbits, params) },
                    "s"    => { circuit.add_gate(S::new(), qbits) },
                    "sdg"  => { circuit.add_gate(Sdg::new(), qbits) },
                    "swap" => { circuit.add_gate(Swap::new(), qbits) },
                    "t"    => { circuit.add_gate(T::new(), qbits) },
                    "tdg"  => { circuit.add_gate(Tdg::new(), qbits) },
                    "u1"   => { add_parametrized_gate!(1, circuit, U1, qbits, params) },
                    "u2"   => { add_parametrized_gate!(2, circuit, U2, qbits, params) },
                    "u3"   => { add_parametrized_gate!(3, circuit, U3, qbits, params) },
                    "v"    => { circuit.add_gate(V::new(), qbits) },
                    "vdg"  => { circuit.add_gate(Vdg::new(), qbits) },
                    "x"    => { circuit.add_gate(X::new(), qbits) },
                    "y"    => { circuit.add_gate(Y::new(), qbits) },
                    "z"    => { circuit.add_gate(Z::new(), qbits) },
                    _ => {
                        Err(crate::error::Error::from(
                            crate::error::ParseError::UnknownGate(String::from(gate_name))
                        ))
                    }
                };

            match res
            {
                Ok(_) => CResult::new(),
                Err(err) => CResult::error(&err.to_string())
            }
        }
        else
        {
            CResult::error("Invalid gate name")
        }
    }
}

#[no_mangle]
pub extern "C" fn circuit_add_conditional_gate(ptr: *mut Circuit,
    control_ptr: *const usize, nr_control: usize, target: u64,
    gate: *const c_char,
    qbits_ptr: *const usize, nr_qbits: usize,
    param_ptr: *const CParameter, nr_params: usize) -> CResult
{
    if ptr.is_null()
    {
        CResult::error("Pointer to circuit is NULL")
    }
    else if control_ptr.is_null()
    {
        CResult::error("Pointer to control bit indices is NULL")
    }
    else if qbits_ptr.is_null()
    {
        CResult::error("Pointer to bit indices is NULL")
    }
    else
    {
        let circuit = unsafe { &mut *ptr };
        let control = unsafe { ::std::slice::from_raw_parts(control_ptr, nr_control) };
        let qbits = unsafe { ::std::slice::from_raw_parts(qbits_ptr, nr_qbits) };
        let params = if param_ptr.is_null()
            {
                &[]
            }
            else
            {
                unsafe { ::std::slice::from_raw_parts(param_ptr, nr_params) }
            };

        if let Ok(gate_name) = unsafe { ::std::ffi::CStr::from_ptr(gate).to_str() }
        {
            let res = match gate_name.to_lowercase().as_str()
            {
                "ch"   => { circuit.add_conditional_gate(control, target, CH::new(), qbits) },
                "cx"   => { circuit.add_conditional_gate(control, target, CX::new(), qbits) },
                "cy"   => { circuit.add_conditional_gate(control, target, CY::new(), qbits) },
                "cz"   => { circuit.add_conditional_gate(control, target, CZ::new(), qbits) },
                "h"    => { circuit.add_conditional_gate(control, target, H::new(), qbits) },
                "i"    => { circuit.add_conditional_gate(control, target, I::new(), qbits) },
                "rx"   => { add_parametrized_gate!(1, circuit, control, target, RX, qbits, params) },
                "ry"   => { add_parametrized_gate!(1, circuit, control, target, RY, qbits, params) },
                "rz"   => { add_parametrized_gate!(1, circuit, control, target, RZ, qbits, params) },
                "s"    => { circuit.add_conditional_gate(control, target, S::new(), qbits) },
                "sdg"  => { circuit.add_conditional_gate(control, target, Sdg::new(), qbits) },
                "swap" => { circuit.add_conditional_gate(control, target, Swap::new(), qbits) },
                "t"    => { circuit.add_conditional_gate(control, target, T::new(), qbits) },
                "tdg"  => { circuit.add_conditional_gate(control, target, Tdg::new(), qbits) },
                "u1"   => { add_parametrized_gate!(1, circuit, control, target, U1, qbits, params) },
                "u2"   => { add_parametrized_gate!(2, circuit, control, target, U2, qbits, params) },
                "u3"   => { add_parametrized_gate!(3, circuit, control, target, U3, qbits, params) },
                "v"    => { circuit.add_conditional_gate(control, target, V::new(), qbits) },
                "vdg"  => { circuit.add_conditional_gate(control, target, Vdg::new(), qbits) },
                "x"    => { circuit.add_conditional_gate(control, target, X::new(), qbits) },
                "y"    => { circuit.add_conditional_gate(control, target, Y::new(), qbits) },
                "z"    => { circuit.add_conditional_gate(control, target, Z::new(), qbits) },
                _ => {
                    Err(crate::error::Error::from(
                        crate::error::ParseError::UnknownGate(String::from(gate_name))
                    ))
                }
            };

            match res
            {
                Ok(_) => CResult::new(),
                Err(err) => CResult::error(&err.to_string())
            }
        }
        else
        {
            CResult::error("Invalid gate name")
        }
    }
}

#[no_mangle]
pub extern "C" fn circuit_reset(ptr: *mut Circuit, qbit: usize) -> CResult
{
    if ptr.is_null()
    {
        CResult::error("Pointer to circuit is NULL")
    }
    else
    {
        let circuit = unsafe { &mut *ptr };
        match circuit.reset(qbit)
        {
            Ok(_) => CResult::new(),
            Err(err) => CResult::error(&err.to_string())
        }
    }
}

#[no_mangle]
pub extern "C" fn circuit_reset_all(ptr: *mut Circuit) -> CResult
{
    if ptr.is_null()
    {
        CResult::error("Pointer to circuit is NULL")
    }
    else
    {
        let circuit = unsafe { &mut *ptr };
        circuit.reset_all();
        CResult::new()
    }
}

#[no_mangle]
pub extern "C" fn circuit_measure(ptr: *mut Circuit, qbit: usize, cbit: usize,
    dir: c_char, collapse: u8) -> CResult
{
    if ptr.is_null()
    {
        CResult::error("Pointer to circuit is NULL")
    }
    else
    {
        let circuit = unsafe { &mut *ptr };
        let basis;
        match dir as u8 as char
        {
            'x'|'X' => { basis = Basis::X; },
            'y'|'Y' => { basis = Basis::Y; },
            'z'|'Z' => { basis = Basis::Z; },
            _       => { return CResult::error(&format!("Invalid measurement basis '{}'", dir)); }
        }

        let res = if collapse != 0 { circuit.measure_basis(qbit, cbit, basis) }
            else { circuit.peek_basis(qbit, cbit, basis) };
        match res
        {
            Ok(_) => CResult::new(),
            Err(err) => CResult::error(&err.to_string())
        }
    }
}

#[no_mangle]
pub extern "C" fn circuit_measure_all(ptr: *mut Circuit,
    cbits_ptr: *const usize, nr_cbits: usize, dir: c_char, collapse: u8) -> CResult
{
    if ptr.is_null()
    {
        CResult::error("Pointer to circuit is NULL")
    }
    else if cbits_ptr.is_null()
    {
        CResult::error("Pointer to measurement bit indices is NULL")
    }
    else
    {
        let cbits = unsafe { ::std::slice::from_raw_parts(cbits_ptr, nr_cbits) };
        let circuit = unsafe { &mut *ptr };
        let basis;
        match dir as u8 as char
        {
            'x'|'X' => { basis = Basis::X; },
            'y'|'Y' => { basis = Basis::Y; },
            'z'|'Z' => { basis = Basis::Z; },
            _       => { return CResult::error(&format!("Invalid measurement basis '{}'", dir)); }
        }

        let res = if collapse != 0 { circuit.measure_all_basis(cbits, basis) }
            else { circuit.peek_all_basis(cbits, basis) };
        match res
        {
            Ok(_) => CResult::new(),
            Err(err) => CResult::error(&err.to_string())
        }
    }
}

#[no_mangle]
pub extern "C" fn circuit_execute(ptr: *mut Circuit, nr_shots: usize) -> CResult
{
    if ptr.is_null()
    {
        CResult::error("Pointer to circuit is NULL")
    }
    else
    {
        let circuit = unsafe { &mut *ptr };
        match circuit.execute(nr_shots)
        {
            Ok(_) => CResult::new(),
            Err(err) => CResult::error(&err.to_string())
        }
    }
}

#[no_mangle]
pub extern "C" fn circuit_reexecute(ptr: *mut Circuit) -> CResult
{
    if ptr.is_null()
    {
        CResult::error("Pointer to circuit is NULL")
    }
    else
    {
        let circuit = unsafe { &mut *ptr };
        match circuit.reexecute()
        {
            Ok(_) => CResult::new(),
            Err(err) => CResult::error(&err.to_string())
        }
    }
}

#[no_mangle]
pub extern "C" fn circuit_histogram(ptr: *const Circuit) -> CResult
{
    if ptr.is_null()
    {
        CResult::error("Pointer to circuit is NULL")
    }
    else
    {
        let circuit = unsafe { &*ptr };
        match circuit.histogram_string()
        {
            Ok(hist) => CResult::histogram(&hist),
            Err(err) => CResult::error(&err.to_string())
        }
    }
}

#[no_mangle]
pub extern "C" fn circuit_latex(ptr: *const Circuit) -> CResult
{
    if ptr.is_null()
    {
        CResult::error("Pointer to circuit is NULL")
    }
    else
    {
        let circuit = unsafe { &*ptr };
        match circuit.latex()
        {
            Ok(ltx) => CResult::string(&ltx),
            Err(err) => CResult::error(&err.to_string())
        }
    }
}


#[no_mangle]
pub extern "C" fn circuit_open_qasm(ptr: *const Circuit) -> CResult
{
    if ptr.is_null()
    {
        CResult::error("Pointer to circuit is NULL")
    }
    else
    {
        let circuit = unsafe { &*ptr };
        match circuit.open_qasm()
        {
            Ok(asm) => CResult::string(&asm),
            Err(err) => CResult::error(&err.to_string())
        }
    }
}

#[no_mangle]
pub extern "C" fn circuit_c_qasm(ptr: *const Circuit) -> CResult
{
    if ptr.is_null()
    {
        CResult::error("Pointer to circuit is NULL")
    }
    else
    {
        let circuit = unsafe { &*ptr };
        match circuit.c_qasm()
        {
            Ok(asm) => CResult::string(&asm),
            Err(err) => CResult::error(&err.to_string())
        }
    }
}
