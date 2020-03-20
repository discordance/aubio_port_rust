use ::libc;
extern "C" {
    #[no_mangle]
    fn calloc(_: libc::c_ulong, _: libc::c_ulong) -> *mut libc::c_void;
    #[no_mangle]
    fn free(_: *mut libc::c_void);
    /* file interface */
    /* strings */
    /* Error reporting */
    /* Logging */
    /* * internal logging function, defined in utils/log.c */
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
pub type aubio_log_level = libc::c_uint;
pub const AUBIO_LOG_LAST_LEVEL: aubio_log_level = 5;
pub const AUBIO_LOG_WRN: aubio_log_level = 4;
pub const AUBIO_LOG_DBG: aubio_log_level = 3;
pub const AUBIO_LOG_MSG: aubio_log_level = 2;
pub const AUBIO_LOG_INF: aubio_log_level = 1;
pub const AUBIO_LOG_ERR: aubio_log_level = 0;
#[derive(Copy, Clone)]
#[repr(C)]
pub struct fmat_t {
    pub length: uint_t,
    pub height: uint_t,
    pub data: *mut *mut smpl_t,
}
#[derive(Copy, Clone)]
#[repr(C)]
pub struct _aubio_sink_t {
    pub sink: *mut libc::c_void,
    pub s_do: aubio_sink_do_t,
    pub s_do_multi: aubio_sink_do_multi_t,
    pub s_preset_samplerate: aubio_sink_preset_samplerate_t,
    pub s_preset_channels: aubio_sink_preset_channels_t,
    pub s_get_samplerate: aubio_sink_get_samplerate_t,
    pub s_get_channels: aubio_sink_get_channels_t,
    pub s_close: aubio_sink_close_t,
    pub s_del: del_aubio_sink_t,
}
pub type del_aubio_sink_t
    =
    Option<unsafe extern "C" fn(_: *mut aubio_sink_t) -> ()>;
/*
  Copyright (C) 2012-2014 Paul Brossier <piem@aubio.org>

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

  Media sink to write blocks of consecutive audio samples to file.

  To read from file, use ::aubio_source_t.

  Depending on how aubio was compiled, the following sinks will be available.

  When creating a new sink using ::new_aubio_sink, the new function of each of
  the compiled-in sinks will be attempted, in the following order, until one of
  them gets successfully created. If all sinks returned NULL, ::new_aubio_sink
  will return NULL.

  \b \p sink_apple_audio : ExtAudioFileRef

  This sink uses CoreAudio [Extended Audio File Services]
  (https://developer.apple.com/library/mac/documentation/MusicAudio/Reference/ExtendedAudioFileServicesReference/Reference/reference.html)
  to write 16-bits encoded WAV files.

  \b \p sink_sndfile : libsndfile

  This sink uses [libsndfile](http://www.mega-nerd.com/libsndfile/) to write
  16-bits encoded WAV files.

  \b \p sink_wavwrite : native WAV write

  A simple sink to write 16-bits PCM RIFF encoded WAV files.

  \example io/test-sink.c

*/
/* * media sink object */
pub type aubio_sink_t = _aubio_sink_t;
pub type aubio_sink_close_t
    =
    Option<unsafe extern "C" fn(_: *mut aubio_sink_t) -> uint_t>;
pub type aubio_sink_get_channels_t
    =
    Option<unsafe extern "C" fn(_: *mut aubio_sink_t) -> uint_t>;
pub type aubio_sink_get_samplerate_t
    =
    Option<unsafe extern "C" fn(_: *mut aubio_sink_t) -> uint_t>;
pub type aubio_sink_preset_channels_t
    =
    Option<unsafe extern "C" fn(_: *mut aubio_sink_t, _: uint_t) -> uint_t>;
pub type aubio_sink_preset_samplerate_t
    =
    Option<unsafe extern "C" fn(_: *mut aubio_sink_t, _: uint_t) -> uint_t>;
pub type aubio_sink_do_multi_t
    =
    Option<unsafe extern "C" fn(_: *mut aubio_sink_t, _: *mut fmat_t,
                                _: uint_t) -> ()>;
/*
  Copyright (C) 2012-2014 Paul Brossier <piem@aubio.org>

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
/* HAVE_SINK_APPLE_AUDIO */
pub type aubio_sink_do_t
    =
    Option<unsafe extern "C" fn(_: *mut aubio_sink_t, _: *mut fvec_t,
                                _: uint_t) -> ()>;
/* *

  create new ::aubio_sink_t

  \param uri the file path or uri to write to
  \param samplerate sample rate to write the file at

  \return newly created ::aubio_sink_t

  Creates a new sink object.

  If samplerate is set to 0, the creation of the file will be delayed until
  both ::aubio_sink_preset_samplerate and ::aubio_sink_preset_channels have
  been called.

*/
#[no_mangle]
pub unsafe extern "C" fn new_aubio_sink(mut uri: *const char_t,
                                        mut samplerate: uint_t)
 -> *mut aubio_sink_t {
    let mut s: *mut aubio_sink_t =
        calloc(::std::mem::size_of::<aubio_sink_t>() as libc::c_ulong,
               1 as libc::c_int as libc::c_ulong) as *mut aubio_sink_t;
    /* HAVE_SINK_APPLE_AUDIO */
    /* HAVE_SNDFILE */
    /* HAVE_WAVWRITE */
    aubio_log(AUBIO_LOG_ERR as libc::c_int,
              b"AUBIO ERROR: sink: failed creating \'%s\' at %dHz (no sink built-in)\n\x00"
                  as *const u8 as *const libc::c_char, uri, samplerate);
    del_aubio_sink(s);
    return 0 as *mut aubio_sink_t;
}
/* *

  write monophonic vector of length hop_size to sink

  \param s sink, created with ::new_aubio_sink
  \param write_data ::fvec_t samples to write to sink
  \param write number of frames to write

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_sink_do(mut s: *mut aubio_sink_t,
                                       mut write_data: *mut fvec_t,
                                       mut write: uint_t) {
    (*s).s_do.expect("non-null function pointer")((*s).sink as
                                                      *mut aubio_sink_t,
                                                  write_data, write);
}
/* *

  write polyphonic vector of length hop_size to sink

  \param s sink, created with ::new_aubio_sink
  \param write_data ::fmat_t samples to write to sink
  \param write number of frames to write

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_sink_do_multi(mut s: *mut aubio_sink_t,
                                             mut write_data: *mut fmat_t,
                                             mut write: uint_t) {
    (*s).s_do_multi.expect("non-null function pointer")((*s).sink as
                                                            *mut aubio_sink_t,
                                                        write_data, write);
}
/* *

  preset sink samplerate

  \param s sink, created with ::new_aubio_sink
  \param samplerate samplerate to preset the sink to, in Hz

  \return 0 on success, 1 on error

  Preset the samplerate of the sink. The file should have been created using a
  samplerate of 0.

  The file will be opened only when both samplerate and channels have been set.

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_sink_preset_samplerate(mut s:
                                                          *mut aubio_sink_t,
                                                      mut samplerate: uint_t)
 -> uint_t {
    return (*s).s_preset_samplerate.expect("non-null function pointer")((*s).sink
                                                                            as
                                                                            *mut aubio_sink_t,
                                                                        samplerate);
}
/* *

  preset sink channels

  \param s sink, created with ::new_aubio_sink
  \param channels number of channels to preset the sink to

  \return 0 on success, 1 on error

  Preset the samplerate of the sink. The file should have been created using a
  samplerate of 0.

  The file will be opened only when both samplerate and channels have been set.

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_sink_preset_channels(mut s: *mut aubio_sink_t,
                                                    mut channels: uint_t)
 -> uint_t {
    return (*s).s_preset_channels.expect("non-null function pointer")((*s).sink
                                                                          as
                                                                          *mut aubio_sink_t,
                                                                      channels);
}
/* *

  get samplerate of sink object

  \param s sink object, created with ::new_aubio_sink
  \return samplerate, in Hz

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_sink_get_samplerate(mut s: *const aubio_sink_t)
 -> uint_t {
    return (*s).s_get_samplerate.expect("non-null function pointer")((*s).sink
                                                                         as
                                                                         *mut aubio_sink_t);
}
/* *

  get channels of sink object

  \param s sink object, created with ::new_aubio_sink
  \return number of channels

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_sink_get_channels(mut s: *const aubio_sink_t)
 -> uint_t {
    return (*s).s_get_channels.expect("non-null function pointer")((*s).sink
                                                                       as
                                                                       *mut aubio_sink_t);
}
/* *

  close sink

  \param s sink object, created with ::new_aubio_sink

  \return 0 on success, non-zero on failure

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_sink_close(mut s: *mut aubio_sink_t)
 -> uint_t {
    return (*s).s_close.expect("non-null function pointer")((*s).sink as
                                                                *mut aubio_sink_t);
}
/* *

  close sink and cleanup memory

  \param s sink object, created with ::new_aubio_sink

*/
#[no_mangle]
pub unsafe extern "C" fn del_aubio_sink(mut s: *mut aubio_sink_t) {
    //AUBIO_ASSERT(s);
    if !s.is_null() && (*s).s_del.is_some() && !(*s).sink.is_null() {
        (*s).s_del.expect("non-null function pointer")((*s).sink as
                                                           *mut aubio_sink_t);
    }
    free(s as *mut libc::c_void);
}
