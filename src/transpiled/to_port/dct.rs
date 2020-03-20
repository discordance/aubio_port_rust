use ::libc;
extern "C" {
    pub type _aubio_dct_plain_t;
    pub type _aubio_dct_ooura_t;
    #[no_mangle]
    fn calloc(_: libc::c_ulong, _: libc::c_ulong) -> *mut libc::c_void;
    #[no_mangle]
    fn free(_: *mut libc::c_void);
    /* * return 1 if a is a power of 2, 0 otherwise */
    #[no_mangle]
    fn aubio_is_power_of_two(a: uint_t) -> uint_t;
    /* Logging */
    /* * internal logging function, defined in utils/log.c */
    #[no_mangle]
    fn aubio_log(level: sint_t, fmt: *const char_t, _: ...) -> uint_t;
    #[no_mangle]
    fn new_aubio_dct_ooura(size: uint_t) -> *mut aubio_dct_ooura_t;
    #[no_mangle]
    fn aubio_dct_ooura_do(s: *mut aubio_dct_ooura_t, input: *const fvec_t,
                          output: *mut fvec_t);
    #[no_mangle]
    fn aubio_dct_ooura_rdo(s: *mut aubio_dct_ooura_t, input: *const fvec_t,
                           output: *mut fvec_t);
    #[no_mangle]
    fn del_aubio_dct_ooura(s: *mut aubio_dct_ooura_t);
    #[no_mangle]
    fn new_aubio_dct_plain(size: uint_t) -> *mut aubio_dct_plain_t;
    #[no_mangle]
    fn aubio_dct_plain_do(s: *mut aubio_dct_plain_t, input: *const fvec_t,
                          output: *mut fvec_t);
    #[no_mangle]
    fn aubio_dct_plain_rdo(s: *mut aubio_dct_plain_t, input: *const fvec_t,
                           output: *mut fvec_t);
    #[no_mangle]
    fn del_aubio_dct_plain(s: *mut aubio_dct_plain_t);
}
/*
  Copyright (C) 2003-2015 Paul Brossier <piem@aubio.org>

  This file is part of aubio.

  aubio is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  aubio is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with aubio.  If not, see <http://www.gnu.org/licenses/>.

*/
/* * \file

  Definition of data types used in aubio

*/
/* * defined to 1 if aubio is compiled in double precision */
/* * short sample format (32 or 64 bits) */
pub type smpl_t = libc::c_float;
/* * print format for sample in double precision */
/* * unsigned integer */
pub type uint_t = libc::c_uint;
/* * signed integer */
pub type sint_t = libc::c_int;
/* * character */
pub type char_t = libc::c_char;
/*
  Copyright (C) 2003-2015 Paul Brossier <piem@aubio.org>

  This file is part of aubio.

  aubio is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  aubio is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with aubio.  If not, see <http://www.gnu.org/licenses/>.

*/
/* * \file

  Vector of real-valued data

  This file specifies the ::fvec_t buffer type, which is used throughout aubio
  to store vector of real-valued ::smpl_t.

  \example test-fvec.c

*/
/* * Buffer for real data

  Vector of real-valued data

  ::fvec_t is is the structure used to store vector of real-valued data, ::smpl_t .

  \code

  uint_t buffer_size = 1024;

  // create a vector of 512 values
  fvec_t * input = new_fvec (buffer_size);

  // set some values of the vector
  input->data[23] = 2.;
  // ..

  // compute the mean of the vector
  mean = fvec_mean(a_vector);

  // destroy the vector
  del_fvec(a_vector);

  \endcode

  See `examples/` and `tests/src` directories for more examples.

 */
#[derive(Copy, Clone)]
#[repr(C)]
pub struct fvec_t {
    pub length: uint_t,
    pub data: *mut smpl_t,
}
/*
  Copyright (C) 2016 Paul Brossier <piem@aubio.org>

  This file is part of aubio.

  aubio is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  aubio is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with aubio.  If not, see <http://www.gnu.org/licenses/>.

*/
/* * \file

  Logging features

  This file specifies ::aubio_log_set_function and
  ::aubio_log_set_level_function, which let you define one or several custom
  logging functions to redirect warnings and errors from aubio to your
  application. The custom function should have the prototype defined in
  ::aubio_log_function_t.

  After a call to ::aubio_log_set_level_function, ::aubio_log_reset can be used
  to reset each logging functions to the default ones.

  \example utils/test-log.c

*/
/* * list of logging levels */
pub type aubio_log_level = libc::c_uint;
/* *< number of valid levels */
/* *< warnings */
pub const AUBIO_LOG_LAST_LEVEL: aubio_log_level = 5;
/* *< debug messages */
pub const AUBIO_LOG_WRN: aubio_log_level = 4;
/* *< general messages */
pub const AUBIO_LOG_DBG: aubio_log_level = 3;
/* *< infos */
pub const AUBIO_LOG_MSG: aubio_log_level = 2;
/* *< critical errors */
pub const AUBIO_LOG_INF: aubio_log_level = 1;
pub const AUBIO_LOG_ERR: aubio_log_level = 0;
#[derive(Copy, Clone)]
#[repr(C)]
pub struct _aubio_dct_t {
    pub dct: *mut libc::c_void,
    pub dct_do: aubio_dct_do_t,
    pub dct_rdo: aubio_dct_rdo_t,
    pub del_dct: del_aubio_dct_t,
}
pub type del_aubio_dct_t
    =
    Option<unsafe extern "C" fn(_: *mut aubio_dct_t) -> ()>;
/*
  Copyright (C) 2017 Paul Brossier <piem@aubio.org>

  This file is part of aubio.

  aubio is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  aubio is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with aubio.  If not, see <http://www.gnu.org/licenses/>.

*/
/* * \file

  Discrete Cosine Transform

  Functions aubio_dct_do() and aubio_dct_rdo() are equivalent to MATLAB/Octave
  dct() and idct() functions, as well as scipy.fftpack.dct(x, norm='ortho') and
  scipy.fftpack.idct(x, norm='ortho')

  \example spectral/test-dct.c

*/
/* * DCT object

  This object computes forward and backward DCT type 2 with orthonormal
  scaling.

*/
pub type aubio_dct_t = _aubio_dct_t;
pub type aubio_dct_rdo_t
    =
    Option<unsafe extern "C" fn(_: *mut aubio_dct_t, _: *const fvec_t,
                                _: *mut fvec_t) -> ()>;
/*
  Copyright (C) 2018 Paul Brossier <piem@aubio.org>

  This file is part of aubio.

  aubio is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  aubio is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with aubio.  If not, see <http://www.gnu.org/licenses/>.

*/
/* * \file

  Discrete Cosine Transform

  Functions aubio_dct_do() and aubio_dct_rdo() are equivalent to MATLAB/Octave
  dct() and idct() functions, as well as scipy.fftpack.dct(x, norm='ortho') and
  scipy.fftpack.idct(x, norm='ortho')

  \example spectral/test-dct.c

*/
// function pointers prototypes
pub type aubio_dct_do_t
    =
    Option<unsafe extern "C" fn(_: *mut aubio_dct_t, _: *const fvec_t,
                                _: *mut fvec_t) -> ()>;
// plain mode
pub type aubio_dct_plain_t = _aubio_dct_plain_t;
pub type aubio_dct_ooura_t = _aubio_dct_ooura_t;
/* * create new DCT computation object

  \param size length of the DCT

*/
#[no_mangle]
pub unsafe extern "C" fn new_aubio_dct(mut size: uint_t) -> *mut aubio_dct_t {
    let mut s: *mut aubio_dct_t =
        calloc(::std::mem::size_of::<aubio_dct_t>() as libc::c_ulong,
               1 as libc::c_int as libc::c_ulong) as *mut aubio_dct_t;
    // ooura support sizes that are power of 2
    if !(aubio_is_power_of_two(size) != 1 as libc::c_int as libc::c_uint ||
             size == 1 as libc::c_int as libc::c_uint) {
        (*s).dct = new_aubio_dct_ooura(size) as *mut libc::c_void;
        if !(*s).dct.is_null() {
            (*s).dct_do =
                ::std::mem::transmute::<Option<unsafe extern "C" fn(_:
                                                                        *mut aubio_dct_ooura_t,
                                                                    _:
                                                                        *const fvec_t,
                                                                    _:
                                                                        *mut fvec_t)
                                                   -> ()>,
                                        aubio_dct_do_t>(Some(aubio_dct_ooura_do
                                                                 as
                                                                 unsafe extern "C" fn(_:
                                                                                          *mut aubio_dct_ooura_t,
                                                                                      _:
                                                                                          *const fvec_t,
                                                                                      _:
                                                                                          *mut fvec_t)
                                                                     -> ()));
            (*s).dct_rdo =
                ::std::mem::transmute::<Option<unsafe extern "C" fn(_:
                                                                        *mut aubio_dct_ooura_t,
                                                                    _:
                                                                        *const fvec_t,
                                                                    _:
                                                                        *mut fvec_t)
                                                   -> ()>,
                                        aubio_dct_rdo_t>(Some(aubio_dct_ooura_rdo
                                                                  as
                                                                  unsafe extern "C" fn(_:
                                                                                           *mut aubio_dct_ooura_t,
                                                                                       _:
                                                                                           *const fvec_t,
                                                                                       _:
                                                                                           *mut fvec_t)
                                                                      -> ()));
            (*s).del_dct =
                ::std::mem::transmute::<Option<unsafe extern "C" fn(_:
                                                                        *mut aubio_dct_ooura_t)
                                                   -> ()>,
                                        del_aubio_dct_t>(Some(del_aubio_dct_ooura
                                                                  as
                                                                  unsafe extern "C" fn(_:
                                                                                           *mut aubio_dct_ooura_t)
                                                                      -> ()));
            return s
        }
        // falling back to plain mode
        aubio_log(AUBIO_LOG_WRN as libc::c_int,
                  b"AUBIO WARNING: dct: no optimised implementation could be created for size %d\n\x00"
                      as *const u8 as *const libc::c_char, size);
    }
    (*s).dct = new_aubio_dct_plain(size) as *mut libc::c_void;
    if !(*s).dct.is_null() {
        (*s).dct_do =
            ::std::mem::transmute::<Option<unsafe extern "C" fn(_:
                                                                    *mut aubio_dct_plain_t,
                                                                _:
                                                                    *const fvec_t,
                                                                _:
                                                                    *mut fvec_t)
                                               -> ()>,
                                    aubio_dct_do_t>(Some(aubio_dct_plain_do as
                                                             unsafe extern "C" fn(_:
                                                                                      *mut aubio_dct_plain_t,
                                                                                  _:
                                                                                      *const fvec_t,
                                                                                  _:
                                                                                      *mut fvec_t)
                                                                 -> ()));
        (*s).dct_rdo =
            ::std::mem::transmute::<Option<unsafe extern "C" fn(_:
                                                                    *mut aubio_dct_plain_t,
                                                                _:
                                                                    *const fvec_t,
                                                                _:
                                                                    *mut fvec_t)
                                               -> ()>,
                                    aubio_dct_rdo_t>(Some(aubio_dct_plain_rdo
                                                              as
                                                              unsafe extern "C" fn(_:
                                                                                       *mut aubio_dct_plain_t,
                                                                                   _:
                                                                                       *const fvec_t,
                                                                                   _:
                                                                                       *mut fvec_t)
                                                                  -> ()));
        (*s).del_dct =
            ::std::mem::transmute::<Option<unsafe extern "C" fn(_:
                                                                    *mut aubio_dct_plain_t)
                                               -> ()>,
                                    del_aubio_dct_t>(Some(del_aubio_dct_plain
                                                              as
                                                              unsafe extern "C" fn(_:
                                                                                       *mut aubio_dct_plain_t)
                                                                  -> ()));
        return s
    } else {
        aubio_log(AUBIO_LOG_ERR as libc::c_int,
                  b"AUBIO ERROR: dct: failed creating with size %d, should be > 0\n\x00"
                      as *const u8 as *const libc::c_char, size);
        del_aubio_dct(s);
        return 0 as *mut aubio_dct_t
    };
}
/* * delete DCT object

  \param s dct object as returned by new_aubio_dct

*/
#[no_mangle]
pub unsafe extern "C" fn del_aubio_dct(mut s: *mut aubio_dct_t) {
    if !(*s).dct.is_null() && (*s).del_dct.is_some() {
        (*s).del_dct.expect("non-null function pointer")((*s).dct as
                                                             *mut aubio_dct_t);
    }
    free(s as *mut libc::c_void);
}
/* * compute forward DCT

  \param s dct object as returned by new_aubio_dct
  \param input input signal
  \param dct_output transformed input array

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_dct_do(mut s: *mut aubio_dct_t,
                                      mut input: *const fvec_t,
                                      mut output: *mut fvec_t) {
    (*s).dct_do.expect("non-null function pointer")((*s).dct as
                                                        *mut aubio_dct_t,
                                                    input, output);
}
/* * compute backward DCT

  \param s dct object as returned by new_aubio_dct
  \param input input signal
  \param idct_output transformed input array

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_dct_rdo(mut s: *mut aubio_dct_t,
                                       mut input: *const fvec_t,
                                       mut output: *mut fvec_t) {
    (*s).dct_rdo.expect("non-null function pointer")((*s).dct as
                                                         *mut aubio_dct_t,
                                                     input, output);
}
