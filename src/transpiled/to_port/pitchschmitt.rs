use ::libc;
extern "C" {
    #[no_mangle]
    fn calloc(_: libc::c_ulong, _: libc::c_ulong) -> *mut libc::c_void;
    #[no_mangle]
    fn free(_: *mut libc::c_void);
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
#[derive(Copy, Clone)]
#[repr(C)]
pub struct _aubio_pitchschmitt_t {
    pub blockSize: uint_t,
    pub rate: uint_t,
    pub schmittBuffer: *mut libc::c_short,
    pub schmittPointer: *mut libc::c_short,
    pub buf: *mut libc::c_short,
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

   Pitch detection using a Schmitt trigger

   This pitch extraction method implements a Schmitt trigger to estimate the
   period of a signal.

   This file was derived from the tuneit project, written by Mario Lang to
   detect the fundamental frequency of a sound.

   See http://delysid.org/tuneit.html

   \example pitch/test-pitchschmitt.c

*/
/* * pitch detection object */
pub type aubio_pitchschmitt_t = _aubio_pitchschmitt_t;
/* * creation of the pitch detection object

  \param buf_size size of the input buffer to analyse

*/
#[no_mangle]
pub unsafe extern "C" fn new_aubio_pitchschmitt(mut size: uint_t)
 -> *mut aubio_pitchschmitt_t {
    let mut p: *mut aubio_pitchschmitt_t =
        calloc(::std::mem::size_of::<aubio_pitchschmitt_t>() as libc::c_ulong,
               1 as libc::c_int as libc::c_ulong) as
            *mut aubio_pitchschmitt_t;
    (*p).blockSize = size;
    (*p).schmittBuffer =
        calloc(((*p).blockSize as
                    libc::c_ulong).wrapping_mul(::std::mem::size_of::<libc::c_short>()
                                                    as libc::c_ulong),
               1 as libc::c_int as libc::c_ulong) as *mut libc::c_short;
    (*p).buf =
        calloc(((*p).blockSize as
                    libc::c_ulong).wrapping_mul(::std::mem::size_of::<libc::c_short>()
                                                    as libc::c_ulong),
               1 as libc::c_int as libc::c_ulong) as *mut libc::c_short;
    (*p).schmittPointer = (*p).schmittBuffer;
    return p;
}
/* * execute pitch detection on an input buffer

  \param p pitch detection object as returned by new_aubio_pitchschmitt
  \param samples_in input signal vector (length as specified at creation time)
  \param cands_out pitch period estimates, in samples

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_pitchschmitt_do(mut p:
                                                   *mut aubio_pitchschmitt_t,
                                               mut input: *const fvec_t,
                                               mut output: *mut fvec_t) {
    let mut j: uint_t = 0;
    j = 0 as libc::c_int as uint_t;
    while j < (*input).length {
        *(*p).buf.offset(j as isize) =
            (*(*input).data.offset(j as isize) as libc::c_double * 32768.0f64)
                as libc::c_short;
        j = j.wrapping_add(1)
    }
    *(*output).data.offset(0 as libc::c_int as isize) =
        aubio_schmittS16LE(p, (*input).length, (*p).buf);
}
/*
  Copyright (C) 2004, 2005  Mario Lang <mlang@delysid.org>
  Copyright (C) 2003-2009 Paul Brossier <piem@aubio.org>

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
#[no_mangle]
pub unsafe extern "C" fn aubio_schmittS16LE(mut p: *mut aubio_pitchschmitt_t,
                                            mut nframes: uint_t,
                                            mut indata: *mut libc::c_short)
 -> smpl_t {
    let mut i: uint_t = 0;
    let mut j: uint_t = 0;
    let mut blockSize: uint_t = (*p).blockSize;
    let mut schmittBuffer: *mut libc::c_short = (*p).schmittBuffer;
    let mut schmittPointer: *mut libc::c_short = (*p).schmittPointer;
    let mut period: smpl_t = 0.0f64 as smpl_t;
    let mut trigfact: smpl_t = 0.6f64 as smpl_t;
    i = 0 as libc::c_int as uint_t;
    while i < nframes {
        let fresh0 = schmittPointer;
        schmittPointer = schmittPointer.offset(1);
        *fresh0 = *indata.offset(i as isize);
        if schmittPointer.wrapping_offset_from(schmittBuffer) as libc::c_long
               >= blockSize as sint_t as libc::c_long {
            let mut endpoint: sint_t = 0;
            let mut startpoint: sint_t = 0;
            let mut t1: sint_t = 0;
            let mut t2: sint_t = 0;
            let mut A1: sint_t = 0;
            let mut A2: sint_t = 0;
            let mut tc: sint_t = 0;
            let mut schmittTriggered: sint_t = 0;
            schmittPointer = schmittBuffer;
            j = 0 as libc::c_int as uint_t;
            A1 = 0 as libc::c_int;
            A2 = 0 as libc::c_int;
            while j < blockSize {
                if *schmittBuffer.offset(j as isize) as libc::c_int >
                       0 as libc::c_int &&
                       A1 < *schmittBuffer.offset(j as isize) as libc::c_int {
                    A1 = *schmittBuffer.offset(j as isize) as sint_t
                }
                if (*schmittBuffer.offset(j as isize) as libc::c_int) <
                       0 as libc::c_int &&
                       A2 <
                           -(*schmittBuffer.offset(j as isize) as libc::c_int)
                   {
                    A2 = -(*schmittBuffer.offset(j as isize) as libc::c_int)
                }
                j = j.wrapping_add(1)
            }
            t1 =
                ((A1 as libc::c_float * trigfact) as libc::c_double + 0.5f64)
                    as sint_t;
            t2 =
                -(((A2 as libc::c_float * trigfact) as libc::c_double +
                       0.5f64) as sint_t);
            startpoint = 0 as libc::c_int;
            j = 1 as libc::c_int as uint_t;
            while j < blockSize &&
                      *schmittBuffer.offset(j as isize) as libc::c_int <= t1 {
                j = j.wrapping_add(1)
            }
            while j < blockSize.wrapping_sub(1 as libc::c_int as libc::c_uint)
                      &&
                      !(*schmittBuffer.offset(j as isize) as libc::c_int >= t2
                            &&
                            (*schmittBuffer.offset(j.wrapping_add(1 as
                                                                      libc::c_int
                                                                      as
                                                                      libc::c_uint)
                                                       as isize) as
                                 libc::c_int) < t2) {
                j = j.wrapping_add(1)
            }
            startpoint = j as sint_t;
            schmittTriggered = 0 as libc::c_int;
            endpoint = startpoint + 1 as libc::c_int;
            j = startpoint as uint_t;
            tc = 0 as libc::c_int;
            while j < blockSize {
                if schmittTriggered == 0 {
                    schmittTriggered =
                        (*schmittBuffer.offset(j as isize) as libc::c_int >=
                             t1) as libc::c_int
                } else if *schmittBuffer.offset(j as isize) as libc::c_int >=
                              t2 &&
                              (*schmittBuffer.offset(j.wrapping_add(1 as
                                                                        libc::c_int
                                                                        as
                                                                        libc::c_uint)
                                                         as isize) as
                                   libc::c_int) < t2 {
                    endpoint = j as sint_t;
                    tc += 1;
                    schmittTriggered = 0 as libc::c_int
                }
                j = j.wrapping_add(1)
            }
            if endpoint > startpoint && tc > 0 as libc::c_int {
                period =
                    (endpoint - startpoint) as smpl_t / tc as libc::c_float
            }
        }
        i = i.wrapping_add(1)
    }
    (*p).schmittBuffer = schmittBuffer;
    (*p).schmittPointer = schmittPointer;
    return period;
}
/* * deletion of the pitch detection object

  \param p pitch detection object as returned by new_aubio_pitchschmitt

*/
#[no_mangle]
pub unsafe extern "C" fn del_aubio_pitchschmitt(mut p:
                                                    *mut aubio_pitchschmitt_t) {
    free((*p).schmittBuffer as *mut libc::c_void);
    free((*p).buf as *mut libc::c_void);
    free(p as *mut libc::c_void);
}
