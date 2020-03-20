extern "C" {
    #[no_mangle]
    fn powf(_: f32, _: f32) -> f32;
    #[no_mangle]
    fn free(_: *mut core::ffi::c_void);
    #[no_mangle]
    fn calloc(_: u64, _: u64) -> *mut core::ffi::c_void;
    /* * fvec_t buffer creation function

  \param length the length of the buffer to create

*/
    #[no_mangle]
    fn new_fvec(length: uint_t) -> *mut fvec_t;
    /* * fvec_t buffer deletion function

  \param s buffer to delete as returned by new_fvec()

*/
    #[no_mangle]
    fn del_fvec(s: *mut fvec_t);
    /* * set all elements to a given value

  \param s vector to modify
  \param val value to set elements to

*/
    #[no_mangle]
    fn fvec_set_all(s: *mut fvec_t, val: smpl_t);
    #[no_mangle]
    fn aubio_log(level: sint_t, fmt: *const char_t, _: ...) -> uint_t;
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
pub type smpl_t = f32;
/* * print format for sample in double precision */
/* * unsigned integer */
pub type uint_t = u32;
/* * signed integer */
pub type sint_t = i32;
/* * character */
pub type char_t = i8;
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
pub type C2RustUnnamed = u32;
pub const AUBIO_FAIL: C2RustUnnamed = 1;
pub const AUBIO_OK: C2RustUnnamed = 0;
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
pub type aubio_log_level = u32;
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
pub struct cvec_t {
    pub length: uint_t,
    pub norm: *mut smpl_t,
    pub phas: *mut smpl_t,
}
// from 1.e-6 to .2
/* * structure to store object state */
#[derive(Copy, Clone)]
#[repr(C)]
pub struct _aubio_spectral_whitening_t {
    pub buf_size: uint_t,
    pub hop_size: uint_t,
    pub samplerate: uint_t,
    pub relax_time: smpl_t,
    pub r_decay: smpl_t,
    pub floor: smpl_t,
    pub peak_values: *mut fvec_t,
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

  Spectral adaptive whitening

  References:

  D. Stowell and M. D. Plumbley. Adaptive whitening for improved real-time
  audio onset detection. In Proceedings of the International Computer Music
  Conference (ICMC), 2007, Copenhagen, Denmark.

  http://www.eecs.qmul.ac.uk/~markp/2007/StowellPlumbley07-icmc.pdf

  S. Böck,, F. Krebs, and M. Schedl. Evaluating the Online Capabilities of
  Onset Detection Methods. In Proceedings of the 13th International Society for
  Music Information Retrieval Conference (ISMIR), 2012, Porto, Portugal.

  http://ismir2012.ismir.net/event/papers/049_ISMIR_2012.pdf
  http://www.cp.jku.at/research/papers/Boeck_etal_ISMIR_2012.pdf

*/
/* * spectral whitening structure */
pub type aubio_spectral_whitening_t = _aubio_spectral_whitening_t;
/* * execute spectral adaptive whitening, in-place

  \param o spectral whitening object as returned by new_aubio_spectral_whitening()
  \param fftgrain input signal spectrum as computed by aubio_pvoc_do() or aubio_fft_do()

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_spectral_whitening_do(o:
                                                         *mut aubio_spectral_whitening_t,
                                                     fftgrain:
                                                         *mut cvec_t) {
    let mut i: uint_t = 0 as i32 as uint_t;
    let length: uint_t =
        if (*fftgrain).length < (*(*o).peak_values).length {
            (*fftgrain).length
        } else { (*(*o).peak_values).length };
    i = 0 as i32 as uint_t;
    while i < length {
        let tmp: smpl_t =
            if (*o).r_decay * *(*(*o).peak_values).data.offset(i as isize) >
                   (*o).floor {
                ((*o).r_decay) * *(*(*o).peak_values).data.offset(i as isize)
            } else { (*o).floor };
        *(*(*o).peak_values).data.offset(i as isize) =
            if *(*fftgrain).norm.offset(i as isize) > tmp {
                *(*fftgrain).norm.offset(i as isize)
            } else { tmp };
        let ref mut fresh0 = *(*fftgrain).norm.offset(i as isize);
        *fresh0 /= *(*(*o).peak_values).data.offset(i as isize);
        i = i.wrapping_add(1)
    };
}
/* * creation of a spectral whitening object

  \param buf_size window size of input grains
  \param hop_size number of samples between two consecutive input grains
  \param samplerate sampling rate of the input signal

*/
#[no_mangle]
pub unsafe extern "C" fn new_aubio_spectral_whitening(buf_size: uint_t,
                                                      hop_size: uint_t,
                                                      samplerate: uint_t)
 -> *mut aubio_spectral_whitening_t {
    let mut o: *mut aubio_spectral_whitening_t =
        calloc(::std::mem::size_of::<aubio_spectral_whitening_t>() as
                   u64, 1 as i32 as u64) as
            *mut aubio_spectral_whitening_t;
    if (buf_size as sint_t) < 1 as i32 {
        aubio_log(AUBIO_LOG_ERR as i32,
                  b"AUBIO ERROR: spectral_whitening: got buffer_size %d, but can not be < 1\n\x00"
                      as *const u8 as *const i8, buf_size);
    } else if (hop_size as sint_t) < 1 as i32 {
        aubio_log(AUBIO_LOG_ERR as i32,
                  b"AUBIO ERROR: spectral_whitening: got hop_size %d, but can not be < 1\n\x00"
                      as *const u8 as *const i8, hop_size);
    } else if (samplerate as sint_t) < 1 as i32 {
        aubio_log(AUBIO_LOG_ERR as i32,
                  b"AUBIO ERROR: spectral_whitening: got samplerate %d, but can not be < 1\n\x00"
                      as *const u8 as *const i8, samplerate);
    } else {
        (*o).peak_values =
            new_fvec(buf_size.wrapping_div(2 as i32 as
                                               u32).wrapping_add(1 as
                                                                              i32
                                                                              as
                                                                              u32));
        (*o).buf_size = buf_size;
        (*o).hop_size = hop_size;
        (*o).samplerate = samplerate;
        (*o).floor = 1.0e-4f64 as smpl_t;
        aubio_spectral_whitening_set_relax_time(o,
                                                250 as i32 as smpl_t);
        aubio_spectral_whitening_reset(o);
        return o
    }
    free(o as *mut core::ffi::c_void);
    return 0 as *mut aubio_spectral_whitening_t;
}
/* * set relaxation time for spectral whitening

  \param o spectral whitening object as returned by new_aubio_spectral_whitening()
  \param relax_time relaxation time in seconds between 20 and 500, defaults 250

  */
#[no_mangle]
pub unsafe extern "C" fn aubio_spectral_whitening_set_relax_time(mut o:
                                                                     *mut aubio_spectral_whitening_t,
                                                                 relax_time:
                                                                     smpl_t)
 -> uint_t {
    (*o).relax_time = relax_time;
    (*o).r_decay =
        powf(0.001f64 as f32,
             (*o).hop_size as f32 / (*o).samplerate as f32
                 / (*o).relax_time);
    return AUBIO_OK as i32 as uint_t;
}
/* * get relaxation time of spectral whitening

  \param o spectral whitening object as returned by new_aubio_spectral_whitening()
  \return relaxation time in seconds

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_spectral_whitening_get_relax_time(o:
                                                                     *mut aubio_spectral_whitening_t)
 -> smpl_t {
    return (*o).relax_time;
}
/* * set floor for spectral whitening

  \param o spectral whitening object as returned by new_aubio_spectral_whitening()
  \param floor value (typically between 1.e-6 and .2, defaults to 1.e-4)

  */
#[no_mangle]
pub unsafe extern "C" fn aubio_spectral_whitening_set_floor(mut o:
                                                                *mut aubio_spectral_whitening_t,
                                                            floor: smpl_t)
 -> uint_t {
    (*o).floor = floor;
    return AUBIO_OK as i32 as uint_t;
}
/* * get floor of spectral whitening

  \param o spectral whitening object as returned by new_aubio_spectral_whitening()
  \return floor value

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_spectral_whitening_get_floor(o:
                                                                *mut aubio_spectral_whitening_t)
 -> smpl_t {
    return (*o).floor;
}
/* * reset spectral whitening object

  \param o spectral whitening object as returned by new_aubio_spectral_whitening()

 */
#[no_mangle]
pub unsafe extern "C" fn aubio_spectral_whitening_reset(o:
                                                            *mut aubio_spectral_whitening_t) {
    /* cover the case n == 0. */
    fvec_set_all((*o).peak_values, (*o).floor);
}
/* * deletion of a spectral whitening

  \param o spectral whitening object as returned by new_aubio_spectral_whitening()

*/
#[no_mangle]
pub unsafe extern "C" fn del_aubio_spectral_whitening(o:
                                                          *mut aubio_spectral_whitening_t) {
    del_fvec((*o).peak_values);
    free(o as *mut core::ffi::c_void);
}
