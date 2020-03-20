use ::libc;
extern "C" {
    pub type _aubio_parameter_t;
    #[no_mangle]
    fn sinf(_: libc::c_float) -> libc::c_float;
    #[no_mangle]
    fn calloc(_: libc::c_ulong, _: libc::c_ulong) -> *mut libc::c_void;
    #[no_mangle]
    fn free(_: *mut libc::c_void);
    #[no_mangle]
    fn floorf(_: libc::c_float) -> libc::c_float;
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
    /* * set all elements to zero

  \param s vector to modify

*/
    #[no_mangle]
    fn fvec_zeros(s: *mut fvec_t);
    /* * clamp the values of a vector within the range [-abs(max), abs(max)]

  \param in vector to clamp
  \param absmax maximum value over which input vector elements should be clamped

*/
    #[no_mangle]
    fn fvec_clamp(in_0: *mut fvec_t, absmax: smpl_t);
    #[no_mangle]
    fn aubio_log(level: sint_t, fmt: *const char_t, _: ...) -> uint_t;
    /* * set all elements to zero

  \param s vector to modify

*/
    #[no_mangle]
    fn fmat_zeros(s: *mut fmat_t);
    /* * create new parameter object

  \param min_value the minimum value of the new parameter
  \param max_value the maximum value of the new parameter
  \param steps the number of steps to interpolate from the old value to the target value

  \return the newly created ::aubio_parameter_t

*/
    #[no_mangle]
    fn new_aubio_parameter(min_value: smpl_t, max_value: smpl_t,
                           steps: uint_t) -> *mut aubio_parameter_t;
    /* * set target value of the parameter

  \param param parameter, created by ::new_aubio_parameter
  \param value new target value

  \return 0 if successful, 1 otherwise

*/
    #[no_mangle]
    fn aubio_parameter_set_target_value(param: *mut aubio_parameter_t,
                                        value: smpl_t) -> uint_t;
    /* * get next parameter

  \param param parameter, created by ::new_aubio_parameter

  \return new interpolated parameter value

*/
    #[no_mangle]
    fn aubio_parameter_get_next_value(param: *mut aubio_parameter_t)
     -> smpl_t;
    /* * get current parameter value, without interpolation

  \param param parameter, created by ::new_aubio_parameter

  \return current value

*/
    #[no_mangle]
    fn aubio_parameter_get_current_value(param: *const aubio_parameter_t)
     -> smpl_t;
    /* * destroy ::aubio_parameter_t object

  \param param parameter, created by ::new_aubio_parameter

*/
    #[no_mangle]
    fn del_aubio_parameter(param: *mut aubio_parameter_t);
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
pub type C2RustUnnamed = libc::c_uint;
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
pub struct fmat_t {
    pub length: uint_t,
    pub height: uint_t,
    pub data: *mut *mut smpl_t,
}
/*
  Copyright (C) 2003-2013 Paul Brossier <piem@aubio.org>

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

  Parameter with linear interpolation

  This object manages a parameter, with minimum and maximum values, and a
  number of steps to compute linear interpolation between two values.

  \example utils/test-parameter.c

*/
/* * parameter object */
pub type aubio_parameter_t = _aubio_parameter_t;
#[derive(Copy, Clone)]
#[repr(C)]
pub struct _aubio_wavetable_t {
    pub samplerate: uint_t,
    pub blocksize: uint_t,
    pub wavetable_length: uint_t,
    pub wavetable: *mut fvec_t,
    pub playing: uint_t,
    pub last_pos: smpl_t,
    pub freq: *mut aubio_parameter_t,
    pub amp: *mut aubio_parameter_t,
}
/*
  Copyright (C) 2003-2013 Paul Brossier <piem@aubio.org>

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

  Wavetable synthesis.

  This file creates a wavetable and plays it at different frequency.

  The `_do` function adds the new samples to the input, and write the result as
  the output.

  \example synth/test-wavetable.c

*/
/* * wavetable object */
pub type aubio_wavetable_t = _aubio_wavetable_t;
/* * create new wavetable object

  \param samplerate the sampling rate of the new wavetable
  \param hop_size the block size of the new wavetable

  \return the newly created aubio_wavetable_t

*/
#[no_mangle]
pub unsafe extern "C" fn new_aubio_wavetable(mut samplerate: uint_t,
                                             mut blocksize: uint_t)
 -> *mut aubio_wavetable_t {
    let mut i: uint_t = 0 as libc::c_int as uint_t;
    let mut s: *mut aubio_wavetable_t =
        calloc(::std::mem::size_of::<aubio_wavetable_t>() as libc::c_ulong,
               1 as libc::c_int as libc::c_ulong) as *mut aubio_wavetable_t;
    if samplerate as sint_t <= 0 as libc::c_int {
        aubio_log(AUBIO_LOG_ERR as libc::c_int,
                  b"AUBIO ERROR: Can not create wavetable with samplerate %d\n\x00"
                      as *const u8 as *const libc::c_char, samplerate);
        free(s as *mut libc::c_void);
        return 0 as *mut aubio_wavetable_t
    } else {
        (*s).samplerate = samplerate;
        (*s).blocksize = blocksize;
        (*s).wavetable_length = 4096 as libc::c_int as uint_t;
        (*s).wavetable =
            new_fvec((*s).wavetable_length.wrapping_add(3 as libc::c_int as
                                                            libc::c_uint));
        i = 0 as libc::c_int as uint_t;
        while i < (*s).wavetable_length {
            *(*(*s).wavetable).data.offset(i as isize) =
                sinf((3.14159265358979323846264338327950288f64 * 2.0f64 *
                          i as libc::c_double /
                          (*s).wavetable_length as smpl_t as libc::c_double)
                         as libc::c_float);
            i = i.wrapping_add(1)
        }
        *(*(*s).wavetable).data.offset((*s).wavetable_length as isize) =
            *(*(*s).wavetable).data.offset(0 as libc::c_int as isize);
        *(*(*s).wavetable).data.offset((*s).wavetable_length.wrapping_add(1 as
                                                                              libc::c_int
                                                                              as
                                                                              libc::c_uint)
                                           as isize) =
            *(*(*s).wavetable).data.offset(1 as libc::c_int as isize);
        *(*(*s).wavetable).data.offset((*s).wavetable_length.wrapping_add(2 as
                                                                              libc::c_int
                                                                              as
                                                                              libc::c_uint)
                                           as isize) =
            *(*(*s).wavetable).data.offset(2 as libc::c_int as isize);
        (*s).playing = 0 as libc::c_int as uint_t;
        (*s).last_pos = 0.0f64 as smpl_t;
        (*s).freq =
            new_aubio_parameter(0.0f64 as smpl_t,
                                ((*s).samplerate as libc::c_double / 2.0f64)
                                    as smpl_t, 10 as libc::c_int as uint_t);
        (*s).amp =
            new_aubio_parameter(0.0f64 as smpl_t, 1.0f64 as smpl_t,
                                100 as libc::c_int as uint_t);
        return s
    };
}
unsafe extern "C" fn interp_2(mut input: *const fvec_t, mut pos: smpl_t)
 -> smpl_t {
    let mut idx: uint_t = floorf(pos) as uint_t;
    let mut frac: smpl_t = pos - idx as smpl_t;
    let mut a: smpl_t = *(*input).data.offset(idx as isize);
    let mut b: smpl_t =
        *(*input).data.offset(idx.wrapping_add(1 as libc::c_int as
                                                   libc::c_uint) as isize);
    return a + frac * (b - a);
}
/* * process wavetable function

  \param o wavetable, created by new_aubio_wavetable()
  \param input input of the wavetable, to be added to the output
  \param output output of the wavetable

This function adds the new samples from the playing wavetable to the output.

If `input` is not NULL and different from `output`, then the samples from `input`
are added to the output.

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_wavetable_do(mut s: *mut aubio_wavetable_t,
                                            mut input: *const fvec_t,
                                            mut output: *mut fvec_t) {
    let mut i: uint_t = 0;
    if (*s).playing != 0 {
        let mut pos: smpl_t = (*s).last_pos;
        i = 0 as libc::c_int as uint_t;
        while i < (*output).length {
            let mut inc: smpl_t = aubio_parameter_get_next_value((*s).freq);
            inc *=
                (*s).wavetable_length as smpl_t / (*s).samplerate as smpl_t;
            pos += inc;
            while pos > (*s).wavetable_length as libc::c_float {
                pos -= (*s).wavetable_length as libc::c_float
            }
            *(*output).data.offset(i as isize) =
                aubio_parameter_get_next_value((*s).amp);
            let ref mut fresh0 = *(*output).data.offset(i as isize);
            *fresh0 *= interp_2((*s).wavetable, pos);
            i = i.wrapping_add(1)
        }
        (*s).last_pos = pos
    } else {
        i = 0 as libc::c_int as uint_t;
        while i < (*output).length {
            aubio_parameter_get_next_value((*s).freq);
            aubio_parameter_get_next_value((*s).amp);
            i = i.wrapping_add(1)
        }
        fvec_zeros(output);
    }
    // add input to output if needed
    if !input.is_null() && input != output as *const fvec_t {
        i = 0 as libc::c_int as uint_t;
        while i < (*output).length {
            let ref mut fresh1 = *(*output).data.offset(i as isize);
            *fresh1 += *(*input).data.offset(i as isize);
            i = i.wrapping_add(1)
        }
        fvec_clamp(output, 1.0f64 as smpl_t);
    };
}
/* * process wavetable function, multiple channels

  \param o wavetable, created by new_aubio_wavetable()
  \param input input of the wavetable, to be added to the output
  \param output output of the wavetable

This function adds the new samples from the playing wavetable to the output.

If `input` is not NULL and different from `output`, then the samples from `input`
are added to the output.

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_wavetable_do_multi(mut s:
                                                      *mut aubio_wavetable_t,
                                                  mut input: *const fmat_t,
                                                  mut output: *mut fmat_t) {
    let mut i: uint_t = 0;
    let mut j: uint_t = 0;
    if (*s).playing != 0 {
        let mut pos: smpl_t = (*s).last_pos;
        j = 0 as libc::c_int as uint_t;
        while j < (*output).length {
            let mut inc: smpl_t = aubio_parameter_get_next_value((*s).freq);
            let mut amp: smpl_t = aubio_parameter_get_next_value((*s).amp);
            inc *=
                (*s).wavetable_length as smpl_t / (*s).samplerate as smpl_t;
            pos += inc;
            while pos > (*s).wavetable_length as libc::c_float {
                pos -= (*s).wavetable_length as libc::c_float
            }
            i = 0 as libc::c_int as uint_t;
            while i < (*output).height {
                *(*(*output).data.offset(i as isize)).offset(j as isize) =
                    amp * interp_2((*s).wavetable, pos);
                i = i.wrapping_add(1)
            }
            j = j.wrapping_add(1)
        }
        (*s).last_pos = pos
    } else {
        j = 0 as libc::c_int as uint_t;
        while j < (*output).length {
            aubio_parameter_get_next_value((*s).freq);
            aubio_parameter_get_next_value((*s).amp);
            j = j.wrapping_add(1)
        }
        fmat_zeros(output);
    }
    // add output to input if needed
    if !input.is_null() && input != output as *const fmat_t {
        i = 0 as libc::c_int as uint_t;
        while i < (*output).height {
            j = 0 as libc::c_int as uint_t;
            while j < (*output).length {
                let ref mut fresh2 =
                    *(*(*output).data.offset(i as isize)).offset(j as isize);
                *fresh2 +=
                    *(*(*input).data.offset(i as isize)).offset(j as isize);
                j = j.wrapping_add(1)
            }
            i = i.wrapping_add(1)
        }
    };
}
/* * get current playing state

  \param o wavetable, created by new_aubio_wavetable()

  \return 0 if not playing, 1 if playing

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_wavetable_get_playing(mut s:
                                                         *const aubio_wavetable_t)
 -> uint_t {
    return (*s).playing;
}
/* * set current playing state

  \param o wavetable, created by new_aubio_wavetable()
  \param playing 0 for not playing, 1 for playing

  \return 0 if successful, 1 otherwise

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_wavetable_set_playing(mut s:
                                                         *mut aubio_wavetable_t,
                                                     mut playing: uint_t)
 -> uint_t {
    (*s).playing =
        if playing == 1 as libc::c_int as libc::c_uint {
            1 as libc::c_int
        } else { 0 as libc::c_int } as uint_t;
    return 0 as libc::c_int as uint_t;
}
/* * play sample from start

  \param o wavetable, created by new_aubio_wavetable()

  \return 0 if successful, 1 otherwise

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_wavetable_play(mut s: *mut aubio_wavetable_t)
 -> uint_t {
    aubio_wavetable_set_amp(s, 0.7f64 as smpl_t);
    return aubio_wavetable_set_playing(s, 1 as libc::c_int as uint_t);
}
/* * stop wavetable

  \param o wavetable, created by new_aubio_wavetable()

  \return 0 if successful, 1 otherwise

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_wavetable_stop(mut s: *mut aubio_wavetable_t)
 -> uint_t {
    //aubio_wavetable_set_freq (s, 0.);
    aubio_wavetable_set_amp(s, 0.0f64 as smpl_t);
    //s->last_pos = 0;
    return aubio_wavetable_set_playing(s, 0 as libc::c_int as uint_t);
}
/* * load source in wavetable

  TODO: This function is not implemented yet. See new_aubio_sampler() instead.

  \param o wavetable, created by new_aubio_wavetable()
  \param uri the uri of the source to load

  \return 0 if successful, non-zero otherwise

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_wavetable_load(mut s: *mut aubio_wavetable_t,
                                              mut uri: *const char_t)
 -> uint_t {
    aubio_log(AUBIO_LOG_ERR as libc::c_int,
              b"AUBIO ERROR: wavetable: load method not implemented yet, see sampler\n\x00"
                  as *const u8 as *const libc::c_char);
    return AUBIO_FAIL as libc::c_int as uint_t;
}
/* * set wavetable frequency

  \param o wavetable, created by new_aubio_wavetable()
  \param freq new frequency value for the wavetable

  \return 0 if successful, 1 otherwise

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_wavetable_set_freq(mut s:
                                                      *mut aubio_wavetable_t,
                                                  mut freq: smpl_t)
 -> uint_t {
    return aubio_parameter_set_target_value((*s).freq, freq);
}
/* * get wavetable frequency

  \param o wavetable, created by new_aubio_wavetable()

  \return current frequency, in Hz

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_wavetable_get_freq(mut s:
                                                      *const aubio_wavetable_t)
 -> smpl_t {
    return aubio_parameter_get_current_value((*s).freq);
}
/* * set wavetable amplitude

  \param o wavetable, created by new_aubio_wavetable()
  \param amp new amplitude value for the wavetable

  \return 0 if successful, 1 otherwise

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_wavetable_set_amp(mut s:
                                                     *mut aubio_wavetable_t,
                                                 mut amp: smpl_t) -> uint_t {
    return aubio_parameter_set_target_value((*s).amp, amp);
}
/* * get wavetable amplitude

  \param o wavetable, created by new_aubio_wavetable()

  \return current amplitude

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_wavetable_get_amp(mut s:
                                                     *const aubio_wavetable_t)
 -> smpl_t {
    return aubio_parameter_get_current_value((*s).amp);
}
/* * destroy aubio_wavetable_t object

  \param o wavetable, created by new_aubio_wavetable()

*/
#[no_mangle]
pub unsafe extern "C" fn del_aubio_wavetable(mut s: *mut aubio_wavetable_t) {
    del_aubio_parameter((*s).freq);
    del_aubio_parameter((*s).amp);
    del_fvec((*s).wavetable);
    free(s as *mut libc::c_void);
}
