use ::libc;
extern "C" {
    pub type _aubio_source_t;
    #[no_mangle]
    fn calloc(_: libc::c_ulong, _: libc::c_ulong) -> *mut libc::c_void;
    #[no_mangle]
    fn free(_: *mut libc::c_void);
    #[no_mangle]
    fn strncpy(_: *mut libc::c_char, _: *const libc::c_char, _: libc::c_ulong)
     -> *mut libc::c_char;
    #[no_mangle]
    fn strnlen(__s1: *const libc::c_char, __n: size_t) -> size_t;
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
    /* Logging */
    /* * internal logging function, defined in utils/log.c */
    #[no_mangle]
    fn aubio_log(level: sint_t, fmt: *const char_t, _: ...) -> uint_t;
    /* * fmat_t buffer creation function

  \param length the length of the matrix to create
  \param height the height of the matrix to create

*/
    #[no_mangle]
    fn new_fmat(height: uint_t, length: uint_t) -> *mut fmat_t;
    /* * fmat_t buffer deletion function

  \param s buffer to delete as returned by new_fmat()

*/
    #[no_mangle]
    fn del_fmat(s: *mut fmat_t);
    /* *

  create new ::aubio_source_t

  \param uri the file path or uri to read from
  \param samplerate sampling rate to view the fie at
  \param hop_size the size of the blocks to read from

  Creates a new source object. If `0` is passed as `samplerate`, the sample
  rate of the original file is used.

  The samplerate of newly created source can be obtained using
  ::aubio_source_get_samplerate.

*/
    #[no_mangle]
    fn new_aubio_source(uri: *const char_t, samplerate: uint_t,
                        hop_size: uint_t) -> *mut aubio_source_t;
    /* *

  read monophonic vector of length hop_size from source object

  \param s source object, created with ::new_aubio_source
  \param read_to ::fvec_t of data to read to
  \param read upon returns, equals to number of frames actually read

  Upon returns, `read` contains the number of frames actually read from the
  source. `hop_size` if enough frames could be read, less otherwise.

*/
    #[no_mangle]
    fn aubio_source_do(s: *mut aubio_source_t, read_to: *mut fvec_t,
                       read: *mut uint_t);
    /* *

  read polyphonic vector of length hop_size from source object

  \param s source object, created with ::new_aubio_source
  \param read_to ::fmat_t of data to read to
  \param[out] read upon returns, equals to number of frames actually read

  Upon returns, `read` contains the number of frames actually read from the
  source. `hop_size` if enough frames could be read, less otherwise.

*/
    #[no_mangle]
    fn aubio_source_do_multi(s: *mut aubio_source_t, read_to: *mut fmat_t,
                             read: *mut uint_t);
    /* *

  seek source object

  \param s source object, created with ::new_aubio_source
  \param pos position to seek to, in frames

  \return 0 if sucessful, non-zero on failure

*/
    #[no_mangle]
    fn aubio_source_seek(s: *mut aubio_source_t, pos: uint_t) -> uint_t;
    /* *

  close source and cleanup memory

  \param s source object, created with ::new_aubio_source

*/
    #[no_mangle]
    fn del_aubio_source(s: *mut aubio_source_t);
}
pub type __darwin_size_t = libc::c_ulong;
pub type size_t = __darwin_size_t;
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
pub struct fmat_t {
    pub length: uint_t,
    pub height: uint_t,
    pub data: *mut *mut smpl_t,
}
/*
  Copyright (C) 2012-2013 Paul Brossier <piem@aubio.org>

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

  Media source to read blocks of consecutive audio samples from file.

  To write to file, use ::aubio_sink_t.

  Depending on how aubio was compiled, the following sources will be available.

  When creating a new source using ::new_aubio_source, the new function of each
  of the compiled-in sources will be used, in the following order, until one of
  them gets successfully created. If all sources returned NULL,
  ::new_aubio_source will return NULL.

  \b \p source_avcodec : libav

  aubio can be optionally compiled with [libav](http://libav.org), which can
  read from a very large number of audio and video formats, including over
  different network protocols such as HTTP.

  \b \p source_apple_audio : ExtAudioFileRef

  On Mac and iOS platforms, aubio should be compiled with CoreAudio [Extended
  Audio File Services]
  (https://developer.apple.com/library/mac/documentation/MusicAudio/Reference/ExtendedAudioFileServicesReference/Reference/reference.html).
  This provides access to most common audio file formats, including compressed
  ones.

  \b \p source_sndfile : libsndfile

  Also optional, aubio can be built against
  [libsndfile](http://www.mega-nerd.com/libsndfile/), which can read [most
  uncompressed formats](http://www.mega-nerd.com/libsndfile/#Features).

  \b \p source_wavread : native WAV reader

  A simple source to read from 16-bits PCM RIFF encoded WAV files.

  \example io/test-source.c

*/
/* * media source object */
pub type aubio_source_t = _aubio_source_t;
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
#[derive(Copy, Clone)]
#[repr(C)]
pub struct _aubio_sampler_t {
    pub samplerate: uint_t,
    pub blocksize: uint_t,
    pub source: *mut aubio_source_t,
    pub source_output: *mut fvec_t,
    pub source_output_multi: *mut fmat_t,
    pub uri: *mut char_t,
    pub playing: uint_t,
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

  Load and play sound files.

  This file loads a sample and gets ready to play it.

  The `_do` function adds the new samples to the input, and write the result as
  the output.

  \example synth/test-sampler.c

*/
/* * sampler object */
pub type aubio_sampler_t = _aubio_sampler_t;
/* * create new sampler object

  \param samplerate the sampling rate of the new sampler
  \param hop_size the block size of the new sampler

  \return the newly created ::aubio_sampler_t

*/
#[no_mangle]
pub unsafe extern "C" fn new_aubio_sampler(mut samplerate: uint_t,
                                           mut blocksize: uint_t)
 -> *mut aubio_sampler_t {
    let mut s: *mut aubio_sampler_t =
        calloc(::std::mem::size_of::<aubio_sampler_t>() as libc::c_ulong,
               1 as libc::c_int as libc::c_ulong) as *mut aubio_sampler_t;
    if (blocksize as sint_t) < 1 as libc::c_int {
        aubio_log(AUBIO_LOG_ERR as libc::c_int,
                  b"AUBIO ERROR: sampler: got blocksize %d, but can not be < 1\n\x00"
                      as *const u8 as *const libc::c_char, blocksize);
        free(s as *mut libc::c_void);
        return 0 as *mut aubio_sampler_t
    } else {
        (*s).samplerate = samplerate;
        (*s).blocksize = blocksize;
        (*s).source_output = new_fvec(blocksize);
        (*s).source_output_multi =
            new_fmat(4 as libc::c_int as uint_t, blocksize);
        (*s).source = 0 as *mut aubio_source_t;
        (*s).playing = 0 as libc::c_int as uint_t;
        return s
    };
}
/* * load source in sampler

  \param o sampler, created by new_aubio_sampler()
  \param uri the uri of the source to load

  \return 0 if successful, non-zero otherwise

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_sampler_load(mut o: *mut aubio_sampler_t,
                                            mut uri: *const char_t)
 -> uint_t {
    if !(*o).source.is_null() { del_aubio_source((*o).source); }
    if !(*o).uri.is_null() { free((*o).uri as *mut libc::c_void); }
    (*o).uri =
        calloc(strnlen(uri,
                       1024 as libc::c_int as
                           size_t).wrapping_mul(::std::mem::size_of::<char_t>()
                                                    as libc::c_ulong),
               1 as libc::c_int as libc::c_ulong) as *mut char_t;
    strncpy((*o).uri, uri, strnlen(uri, 1024 as libc::c_int as size_t));
    (*o).source = new_aubio_source(uri, (*o).samplerate, (*o).blocksize);
    if !(*o).source.is_null() { return 0 as libc::c_int as uint_t }
    aubio_log(AUBIO_LOG_ERR as libc::c_int,
              b"AUBIO ERROR: sampler: failed loading %s\x00" as *const u8 as
                  *const libc::c_char, uri);
    return 1 as libc::c_int as uint_t;
}
/* * process sampler function

  \param o sampler, created by new_aubio_sampler()
  \param input input of the sampler, to be added to the output
  \param output output of the sampler

This function adds the new samples from the playing source to the output.

If `input` is not NULL and different from `output`, then the samples from `input`
are added to the output.

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_sampler_do(mut o: *mut aubio_sampler_t,
                                          mut input: *const fvec_t,
                                          mut output: *mut fvec_t) {
    let mut read: uint_t = 0 as libc::c_int as uint_t;
    let mut i: uint_t = 0;
    if (*o).playing != 0 {
        aubio_source_do((*o).source, (*o).source_output, &mut read);
        i = 0 as libc::c_int as uint_t;
        while i < (*output).length {
            let ref mut fresh0 = *(*output).data.offset(i as isize);
            *fresh0 += *(*(*o).source_output).data.offset(i as isize);
            i = i.wrapping_add(1)
        }
        if read < (*o).blocksize { (*o).playing = 0 as libc::c_int as uint_t }
    }
    if !input.is_null() && input != output as *const fvec_t {
        i = 0 as libc::c_int as uint_t;
        while i < (*output).length {
            let ref mut fresh1 = *(*output).data.offset(i as isize);
            *fresh1 += *(*input).data.offset(i as isize);
            i = i.wrapping_add(1)
        }
    };
}
/* * process sampler function, multiple channels

  \param o sampler, created by new_aubio_sampler()
  \param input input of the sampler, to be added to the output
  \param output output of the sampler

This function adds the new samples from the playing source to the output.

If `input` is not NULL and different from `output`, then the samples from `input`
are added to the output.

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_sampler_do_multi(mut o: *mut aubio_sampler_t,
                                                mut input: *const fmat_t,
                                                mut output: *mut fmat_t) {
    let mut read: uint_t = 0 as libc::c_int as uint_t;
    let mut i: uint_t = 0;
    let mut j: uint_t = 0;
    if (*o).playing != 0 {
        aubio_source_do_multi((*o).source, (*o).source_output_multi,
                              &mut read);
        i = 0 as libc::c_int as uint_t;
        while i < (*output).height {
            j = 0 as libc::c_int as uint_t;
            while j < (*output).length {
                let ref mut fresh2 =
                    *(*(*output).data.offset(i as isize)).offset(j as isize);
                *fresh2 +=
                    *(*(*(*o).source_output_multi).data.offset(i as
                                                                   isize)).offset(j
                                                                                      as
                                                                                      isize);
                j = j.wrapping_add(1)
            }
            i = i.wrapping_add(1)
        }
        if read < (*o).blocksize { (*o).playing = 0 as libc::c_int as uint_t }
    }
    if !input.is_null() && input != output as *const fmat_t {
        i = 0 as libc::c_int as uint_t;
        while i < (*output).height {
            j = 0 as libc::c_int as uint_t;
            while j < (*output).length {
                let ref mut fresh3 =
                    *(*(*output).data.offset(i as isize)).offset(j as isize);
                *fresh3 +=
                    *(*(*input).data.offset(i as isize)).offset(j as isize);
                j = j.wrapping_add(1)
            }
            i = i.wrapping_add(1)
        }
    };
}
/* * get current playing state

  \param o sampler, created by new_aubio_sampler()

  \return 0 if not playing, 1 if playing

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_sampler_get_playing(mut o:
                                                       *const aubio_sampler_t)
 -> uint_t {
    return (*o).playing;
}
/* * set current playing state

  \param o sampler, created by new_aubio_sampler()
  \param playing 0 for not playing, 1 for playing

  \return 0 if successful, 1 otherwise

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_sampler_set_playing(mut o:
                                                       *mut aubio_sampler_t,
                                                   mut playing: uint_t)
 -> uint_t {
    (*o).playing =
        if playing == 1 as libc::c_int as libc::c_uint {
            1 as libc::c_int
        } else { 0 as libc::c_int } as uint_t;
    return 0 as libc::c_int as uint_t;
}
/* * play sample from start

  \param o sampler, created by new_aubio_sampler()

  \return 0 if successful, 1 otherwise

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_sampler_play(mut o: *mut aubio_sampler_t)
 -> uint_t {
    aubio_source_seek((*o).source, 0 as libc::c_int as uint_t);
    return aubio_sampler_set_playing(o, 1 as libc::c_int as uint_t);
}
/* * stop sample

  \param o sampler, created by new_aubio_sampler()

  \return 0 if successful, 1 otherwise

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_sampler_stop(mut o: *mut aubio_sampler_t)
 -> uint_t {
    return aubio_sampler_set_playing(o, 0 as libc::c_int as uint_t);
}
/* * destroy ::aubio_sampler_t object

  \param o sampler, created by new_aubio_sampler()

*/
#[no_mangle]
pub unsafe extern "C" fn del_aubio_sampler(mut o: *mut aubio_sampler_t) {
    if !(*o).source.is_null() { del_aubio_source((*o).source); }
    if !(*o).uri.is_null() { free((*o).uri as *mut libc::c_void); }
    del_fvec((*o).source_output);
    del_fmat((*o).source_output_multi);
    free(o as *mut libc::c_void);
}
