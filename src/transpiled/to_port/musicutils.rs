use ::libc;
extern "C" {
    #[no_mangle]
    fn expf(_: libc::c_float) -> libc::c_float;
    #[no_mangle]
    fn logf(_: libc::c_float) -> libc::c_float;
    #[no_mangle]
    fn powf(_: libc::c_float, _: libc::c_float) -> libc::c_float;
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
pub const AUBIO_LOG_WRN: aubio_log_level = 4;
pub type aubio_log_level = libc::c_uint;
pub const AUBIO_LOG_LAST_LEVEL: aubio_log_level = 5;
pub const AUBIO_LOG_DBG: aubio_log_level = 3;
pub const AUBIO_LOG_MSG: aubio_log_level = 2;
pub const AUBIO_LOG_INF: aubio_log_level = 1;
pub const AUBIO_LOG_ERR: aubio_log_level = 0;
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
#[no_mangle]
pub unsafe extern "C" fn aubio_hztomel(mut freq: smpl_t) -> smpl_t {
    let lin_space: smpl_t = (3.0f64 / 200.0f64) as smpl_t;
    let split_hz: smpl_t = 1000.0f64 as smpl_t;
    let split_mel: smpl_t = split_hz * lin_space;
    let log_space: smpl_t =
        (27.0f64 /
             logf((6400 as libc::c_int as libc::c_double / 1000.0f64) as
                      libc::c_float) as libc::c_double) as smpl_t;
    if freq < 0 as libc::c_int as libc::c_float {
        aubio_log(AUBIO_LOG_WRN as libc::c_int,
                  b"AUBIO WARNING: hztomel: input frequency should be >= 0\n\x00"
                      as *const u8 as *const libc::c_char);
        return 0 as libc::c_int as smpl_t
    }
    if freq < split_hz {
        return freq * lin_space
    } else { return split_mel + log_space * logf(freq / split_hz) };
}
#[no_mangle]
pub unsafe extern "C" fn aubio_meltohz(mut mel: smpl_t) -> smpl_t {
    let lin_space: smpl_t = (200.0f64 / 3.0f64) as smpl_t;
    let split_hz: smpl_t = 1000.0f64 as smpl_t;
    let split_mel: smpl_t = split_hz / lin_space;
    let logSpacing: smpl_t =
        powf((6400 as libc::c_int as libc::c_double / 1000.0f64) as
                 libc::c_float,
             (1 as libc::c_int as libc::c_double / 27.0f64) as libc::c_float);
    if mel < 0 as libc::c_int as libc::c_float {
        aubio_log(AUBIO_LOG_WRN as libc::c_int,
                  b"AUBIO WARNING: meltohz: input mel should be >= 0\n\x00" as
                      *const u8 as *const libc::c_char);
        return 0 as libc::c_int as smpl_t
    }
    if mel < split_mel {
        return lin_space * mel
    } else { return split_hz * powf(logSpacing, mel - split_mel) };
}
#[no_mangle]
pub unsafe extern "C" fn aubio_hztomel_htk(mut freq: smpl_t) -> smpl_t {
    let split_hz: smpl_t = 700.0f64 as smpl_t;
    let log_space: smpl_t = 1127.0f64 as smpl_t;
    if freq < 0 as libc::c_int as libc::c_float {
        aubio_log(AUBIO_LOG_WRN as libc::c_int,
                  b"AUBIO WARNING: hztomel_htk: input frequency should be >= 0\n\x00"
                      as *const u8 as *const libc::c_char);
        return 0 as libc::c_int as smpl_t
    }
    return log_space *
               logf(1 as libc::c_int as libc::c_float + freq / split_hz);
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
 *  various functions useful in audio signal processing
 */
/* * create window

  \param window_type type of the window to create
  \param size length of the window to create (see fvec_set_window())

*/
/* * set elements of a vector to window coefficients

  \param window exsting ::fvec_t to use
  \param window_type type of the window to create

  List of available window types: "rectangle", "hamming", "hanning",
  "hanningz", "blackman", "blackman_harris", "gaussian", "welch", "parzen",
  "default".

  "default" is equivalent to "hanningz".

  References:

    - <a href="http://en.wikipedia.org/wiki/Window_function">Window
function</a> on Wikipedia
    - Amalia de Götzen, Nicolas Bernardini, and Daniel Arfib. Traditional (?)
implementations of a phase vocoder: the tricks of the trade. In Proceedings of
the International Conference on Digital Audio Effects (DAFx-00), pages 37–44,
Uni- versity of Verona, Italy, 2000.
  (<a href="http://www.cs.princeton.edu/courses/archive/spr09/cos325/Bernardini.pdf">
  pdf</a>)

 */
/* * compute the principal argument

  This function maps the input phase to its corresponding value wrapped in the
range \f$ [-\pi, \pi] \f$.

  \param phase unwrapped phase to map to the unit circle

  \return equivalent phase wrapped to the unit circle

*/
/* * convert frequency bin to midi value */
/* * convert midi value to frequency bin */
/* * convert frequency bin to frequency (Hz) */
/* * convert frequency (Hz) to frequency bin */
/* * convert frequency (Hz) to mel

  \param freq input frequency, in Hz

  \return output mel

  Converts a scalar from the frequency domain to the mel scale using Slaney
  Auditory Toolbox's implementation:

  If \f$ f < 1000 \f$, \f$ m = 3 f / 200 \f$.

  If \f$ f >= 1000 \f$, \f$ m = 1000 + 27 \frac{{ln}(f) - ln(1000))}
  {{ln}(6400) - ln(1000)}
  \f$

  See also
  --------

  aubio_meltohz(), aubio_hztomel_htk().

*/
/* * convert mel to frequency (Hz)

  \param mel input mel

  \return output frequency, in Hz

  Converts a scalar from the mel scale to the frequency domain using Slaney
  Auditory Toolbox's implementation:

  If \f$ f < 1000 \f$, \f$ f = 200 m/3 \f$.

  If \f$ f \geq 1000 \f$, \f$ f = 1000 + \left(\frac{6400}{1000}\right)
  ^{\frac{m - 1000}{27}} \f$

  See also
  --------

  aubio_hztomel(), aubio_meltohz_htk().

  References
  ----------

  Malcolm Slaney, *Auditory Toolbox Version 2, Technical Report #1998-010*
  https://engineering.purdue.edu/~malcolm/interval/1998-010/

*/
/* * convert frequency (Hz) to mel

  \param freq input frequency, in Hz

  \return output mel

  Converts a scalar from the frequency domain to the mel scale, using the
  equation defined by O'Shaughnessy, as implemented in the HTK speech
  recognition toolkit:

  \f$ m = 1127 + ln(1 + \frac{f}{700}) \f$

  See also
  --------

  aubio_meltohz_htk(), aubio_hztomel().

  References
  ----------

  Douglas O'Shaughnessy (1987). *Speech communication: human and machine*.
  Addison-Wesley. p. 150. ISBN 978-0-201-16520-3.

  HTK Speech Recognition Toolkit: http://htk.eng.cam.ac.uk/

 */
/* * convert mel to frequency (Hz)

  \param mel input mel

  \return output frequency, in Hz

  Converts a scalar from the mel scale to the frequency domain, using the
  equation defined by O'Shaughnessy, as implemented in the HTK speech
  recognition toolkit:

  \f$ f = 700 * {e}^\left(\frac{f}{1127} - 1\right) \f$

  See also
  --------

  aubio_hztomel_htk(), aubio_meltohz().

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_meltohz_htk(mut mel: smpl_t) -> smpl_t {
    let split_hz: smpl_t = 700.0f64 as smpl_t;
    let log_space: smpl_t = (1.0f64 / 1127.0f64) as smpl_t;
    if mel < 0 as libc::c_int as libc::c_float {
        aubio_log(AUBIO_LOG_WRN as libc::c_int,
                  b"AUBIO WARNING: meltohz_htk: input frequency should be >= 0\n\x00"
                      as *const u8 as *const libc::c_char);
        return 0 as libc::c_int as smpl_t
    }
    return (split_hz as libc::c_double *
                (expf(mel * log_space) as libc::c_double - 1.0f64)) as smpl_t;
}
