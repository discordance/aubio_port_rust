extern "C" {
    #[no_mangle]
    static mut __stdoutp: *mut FILE;
    #[no_mangle]
    static mut __stderrp: *mut FILE;
    #[no_mangle]
    fn fprintf(_: *mut FILE, _: *const i8, _: ...) -> i32;
    #[no_mangle]
    fn vsnprintf(_: *mut i8, _: u64,
                 _: *const i8, _: ::std::ffi::VaList)
     -> i32;
}
pub type __builtin_va_list = [__va_list_tag; 1];
#[derive(Copy, Clone)]
#[repr(C)]
pub struct __va_list_tag {
    pub gp_offset: u32,
    pub fp_offset: u32,
    pub overflow_arg_area: *mut core::ffi::c_void,
    pub reg_save_area: *mut core::ffi::c_void,
}
pub type __int64_t = i64;
pub type __darwin_va_list = __builtin_va_list;
pub type __darwin_off_t = __int64_t;
pub type va_list = __darwin_va_list;
pub type fpos_t = __darwin_off_t;
#[derive(Copy, Clone)]
#[repr(C)]
pub struct __sbuf {
    pub _base: *mut u8,
    pub _size: i32,
}
#[derive(Copy, Clone)]
#[repr(C)]
pub struct __sFILE {
    pub _p: *mut u8,
    pub _r: i32,
    pub _w: i32,
    pub _flags: i16,
    pub _file: i16,
    pub _bf: __sbuf,
    pub _lbfsize: i32,
    pub _cookie: *mut core::ffi::c_void,
    pub _close: Option<unsafe extern "C" fn(_: *mut core::ffi::c_void)
                           -> i32>,
    pub _read: Option<unsafe extern "C" fn(_: *mut core::ffi::c_void,
                                           _: *mut i8,
                                           _: i32) -> i32>,
    pub _seek: Option<unsafe extern "C" fn(_: *mut core::ffi::c_void, _: fpos_t,
                                           _: i32) -> fpos_t>,
    pub _write: Option<unsafe extern "C" fn(_: *mut core::ffi::c_void,
                                            _: *const i8,
                                            _: i32) -> i32>,
    pub _ub: __sbuf,
    pub _ur: i32,
    pub _ubuf: [u8; 3],
    pub _nbuf: [u8; 1],
    pub _lb: __sbuf,
    pub _blksize: i32,
    pub _offset: fpos_t,
}
pub type FILE = __sFILE;
/* * print format for sample in double precision */
/* * unsigned integer */
pub type uint_t = u32;
/* * signed integer */
pub type sint_t = i32;
/* * character */
pub type char_t = i8;
pub type C2RustUnnamed = u32;
pub const AUBIO_FAIL: C2RustUnnamed = 1;
pub const AUBIO_OK: C2RustUnnamed = 0;
pub type aubio_log_level = u32;
pub const AUBIO_LOG_LAST_LEVEL: aubio_log_level = 5;
pub const AUBIO_LOG_WRN: aubio_log_level = 4;
pub const AUBIO_LOG_DBG: aubio_log_level = 3;
pub const AUBIO_LOG_MSG: aubio_log_level = 2;
pub const AUBIO_LOG_INF: aubio_log_level = 1;
pub const AUBIO_LOG_ERR: aubio_log_level = 0;
pub type aubio_log_function_t
    =
    Option<unsafe extern "C" fn(_: sint_t, _: *const char_t,
                                _: *mut core::ffi::c_void) -> ()>;
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
/* * array of pointers to logging functions, one per level */
static mut aubio_log_function: [aubio_log_function_t; 5] = [None; 5];
/* * array of pointers to closure passed to logging functions, one per level */
static mut aubio_log_user_data: [*mut core::ffi::c_void; 5] =
    [0 as *const core::ffi::c_void as *mut core::ffi::c_void; 5];
/* * buffer for logging messages */
static mut aubio_log_buffer: [i8; 512] = [0; 512];
/* * private function used by default by logging functions */
#[no_mangle]
pub unsafe extern "C" fn aubio_default_log(level: sint_t,
                                           message: *const char_t,
                                           _data: *mut core::ffi::c_void) {
    let mut out: *mut FILE = 0 as *mut FILE;
    out = __stdoutp;
    if level == AUBIO_LOG_ERR as i32 ||
           level == AUBIO_LOG_DBG as i32 ||
           level == AUBIO_LOG_WRN as i32 {
        out = __stderrp
    }
    fprintf(out, b"%s\x00" as *const u8 as *const i8, message);
    //fflush(out);
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
/* * @file
 * Private include file
 *
 * This file is for inclusion from _within_ the library only.
 */
/* ********************
 *
 * External includes
 *
 */
/* must be included before fftw3.h */
// for CHAR_BIT, in C99 standard
// --enable-blas=true
/* HAVE_ACCELERATE */
/* HAVE_BLAS */
/* ***
 *
 * SYSTEM INTERFACE
 *
 */
/* Memory management */
/* file interface */
/* strings */
/* Error reporting */
/* Logging */
/* * internal logging function, defined in utils/log.c */
#[no_mangle]
pub unsafe extern "C" fn aubio_log(level: sint_t, fmt: *const char_t,
                                   args: ...) -> uint_t {
    let mut fun: aubio_log_function_t = None;
    let mut args_0: ::std::ffi::VaListImpl;
    args_0 = args.clone();
    vsnprintf(aubio_log_buffer.as_mut_ptr(),
              ::std::mem::size_of::<[i8; 512]>() as u64,
              fmt, args_0.as_va_list());
    if level >= 0 as i32 &&
           level < AUBIO_LOG_LAST_LEVEL as i32 {
        fun = aubio_log_function[level as usize];
        if fun.is_some() {
            Some(fun.expect("non-null function pointer")).expect("non-null function pointer")(level,
                                                                                              aubio_log_buffer.as_mut_ptr(),
                                                                                              aubio_log_user_data[level
                                                                                                                      as
                                                                                                                      usize]);
        } else {
            aubio_default_log(level, aubio_log_buffer.as_mut_ptr(),
                              0 as *mut core::ffi::c_void);
        }
    }
    return AUBIO_FAIL as i32 as uint_t;
}
/* * Reset all logging functions to the default one

 After calling this function, the default logging function will be used to
 print error, warning, normal, and debug messages to `stdout` or `stderr`.

 */
#[no_mangle]
pub unsafe extern "C" fn aubio_log_reset() {
    let mut i: uint_t = 0 as i32 as uint_t;
    i = 0 as i32 as uint_t;
    while i < AUBIO_LOG_LAST_LEVEL as i32 as u32 {
        aubio_log_set_level_function(i as sint_t,
                                     Some(aubio_default_log as
                                              unsafe extern "C" fn(_: sint_t,
                                                                   _:
                                                                       *const char_t,
                                                                   _:
                                                                       *mut core::ffi::c_void)
                                                  -> ()),
                                     0 as *mut core::ffi::c_void);
        i = i.wrapping_add(1)
    };
}
/* * Set logging function for a given level

  \param level the level for which to set the logging function
  \param fun the function to be used to log, of type ::aubio_log_function_t
  \param data optional closure to be passed to the function (can be NULL if
  nothing to pass)

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_log_set_level_function(level: sint_t,
                                                      fun:
                                                          aubio_log_function_t,
                                                      data:
                                                          *mut core::ffi::c_void)
 -> aubio_log_function_t {
    let mut old: aubio_log_function_t = None;
    if level >= 0 as i32 &&
           level < AUBIO_LOG_LAST_LEVEL as i32 {
        old = aubio_log_function[level as usize];
        aubio_log_function[level as usize] = fun;
        aubio_log_user_data[level as usize] = data
    }
    return old;
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
/* *< critical errors */
/* *< infos */
/* *< general messages */
/* *< debug messages */
/* *< warnings */
/* *< number of valid levels */
/* * Logging function prototype, to be passed to ::aubio_log_set_function

  \param level log level
  \param message text to log
  \param data optional closure used by the callback

  See @ref utils/test-log.c for an example of logging function.

 */
/* * Set logging function for all levels

  \param fun the function to be used to log, of type ::aubio_log_function_t
  \param data optional closure to be passed to the function (can be NULL if
  nothing to pass)

 */
#[no_mangle]
pub unsafe extern "C" fn aubio_log_set_function(fun: aubio_log_function_t,
                                                data: *mut core::ffi::c_void) {
    let mut i: uint_t = 0 as i32 as uint_t;
    i = 0 as i32 as uint_t;
    while i < AUBIO_LOG_LAST_LEVEL as i32 as u32 {
        aubio_log_set_level_function(i as sint_t, fun, data);
        i = i.wrapping_add(1)
    };
}
