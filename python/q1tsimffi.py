import cffi

RESULT_ERROR = 0
RESULT_EMPTY = 1
RESULT_STRING = 2
RESULT_HISTOGRAM = 3

ffi = cffi.FFI()
ffi.cdef("""
    typedef struct circuit circuit_t;
    typedef struct
    {
        const char* key;
        size_t count;
    } histelem_t;
    typedef struct
    {
        double value;
        double* value_ptr;
    } parameter_t;
    typedef struct
    {
        void* data;
        size_t length;
        size_t size;
        uint32_t restype;
    } result_t;

    void result_free(result_t res);

    circuit_t *circuit_new(size_t nr_qbits, size_t nr_cbits);
    void circuit_free(circuit_t *ptr);
    size_t circuit_nr_qbits(const circuit_t *ptr);
    size_t circuit_nr_cbits(const circuit_t *ptr);
    result_t circuit_add_gate(circuit_t *ptr, const char* gate,
        const size_t* qbits, size_t nr_qubits,
        const parameter_t* param_ptr, size_t nr_params);
    result_t circuit_add_conditional_gate(circuit_t *ptr,
        const size_t* control_ptr, size_t nr_control, uint64_t target,
        const char* gate,
        const size_t* qbits_ptr, size_t nr_qbits,
        const parameter_t* param_ptr, size_t nr_params);
    result_t circuit_measure(circuit_t *ptr, size_t qbit, size_t cbit, char dir,
        uint8_t collapse);
    result_t circuit_measure_all(circuit_t *ptr, const size_t* cbits, size_t nr_cbits,
        char dir, uint8_t collapse);
    result_t circuit_reset(circuit_t *ptr, size_t qbit);
    result_t circuit_reset_all(circuit_t *ptr);
    result_t circuit_execute(circuit_t *ptr, size_t nr_shots);
    result_t circuit_reexecute(circuit_t *ptr);
    result_t circuit_histogram(const circuit_t *ptr);
    result_t circuit_latex(const circuit_t *ptr);
    result_t circuit_open_qasm(const circuit_t *ptr);
    result_t circuit_c_qasm(const circuit_t *ptr);
""")
qsim_obj = ffi.dlopen("../target/debug/libq1tsim.so")

def unpack_result(res):
    restype = res.restype
    try:
        if restype == RESULT_ERROR:
            msg = ffi.string(ffi.cast("const char*", res.data)).decode('utf-8')
            raise Exception(msg)
        if res.restype == RESULT_EMPTY:
            return None
        elif res.restype == RESULT_STRING:
            msg = ffi.string(ffi.cast("const char*", res.data)).decode('utf-8')
            return msg
        elif res.restype == RESULT_HISTOGRAM:
            elems = ffi.cast("const histelem_t*", res.data)
            hist = dict((ffi.string(elems[i].key).decode('utf-8'), elems[i].count)
                        for i in range(res.length))
            return hist
        else:
            raise Exception("Unknown data type code: {}".format(res.restype))
    finally:
        qsim_obj.result_free(res)

def make_parameters(values):
    res = ffi.new('parameter_t[]', len(values))
    for (i, val) in enumerate(values):
        if type(val) == float:
            res[i].value = val
            res[i].value_ptr = ffi.NULL
        else:
            res[i].value = 0.0
            res[i].value_ptr = val.ptr
    return res

def q1tsim():
    return qsim_obj
