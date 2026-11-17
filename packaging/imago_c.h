/****************************************************************************
 * Copyright (C) from 2009 to Present EPAM Systems.
 *
 * This file is part of Imago toolkit.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 ***************************************************************************/

/*
 * Self-contained consumer header for the Apache-2.0 Imago build vendored at
 * imago/lib/imago.{dll,so,dylib} (built fresh from github.com/epam/Imago, not
 * the upstream api/c/src/imago_c.h, which requires Imago's internal build
 * macros such as CEXPORT/qword/DLLEXPORT to already be defined by the
 * includer). Declares exactly the functions ImagoService.cpp calls, matching
 * the real exported symbol names (verified via objdump against the built
 * Windows DLL; the C API and its exported names are identical on Linux/macOS
 * since all three platforms build from the same upstream source).
 *
 * Tracked here (unlike the imago/ directory itself, which is gitignored like
 * indigo/ -- see .gitignore) because, unlike Indigo's own public headers,
 * nothing upstream ships a header a consumer can use directly: every
 * platform's build_imago.sh (or the manual Windows steps in README.md) just
 * copies this file into imago/include/imago_c.h.
 */

#pragma once
#ifndef _imago_c_h
#define _imago_c_h

#ifdef _WIN32
#define IMAGO_DLLIMPORT __declspec(dllimport)
#else
#define IMAGO_DLLIMPORT
#endif

/* Must match core/src/comdef.h's own qword exactly - MSVC's _int64 keyword is
 * meaningless to GCC even under MinGW, where _WIN32 is defined too. */
#if defined(_WIN32) && !defined(__GNUC__)
#define imago_qword unsigned _int64
#else
#define imago_qword unsigned long long
#endif

#ifndef __cplusplus
#define IMAGO_CEXPORT IMAGO_DLLIMPORT
#else
#define IMAGO_CEXPORT extern "C" IMAGO_DLLIMPORT
#endif

/* All functions return 0 on error, and a nonzero value on success. */

IMAGO_CEXPORT const char *imagoGetVersion();
IMAGO_CEXPORT const char *imagoGetLastError();

IMAGO_CEXPORT imago_qword imagoAllocSessionId();
IMAGO_CEXPORT void imagoSetSessionId(imago_qword id);
IMAGO_CEXPORT void imagoReleaseSessionId(imago_qword id);

IMAGO_CEXPORT int imagoLoadImageFromFile(const char *fileName);
IMAGO_CEXPORT int imagoLoadImageFromBuffer(const char *buf, const int buf_size);
IMAGO_CEXPORT int imagoFilterImage();

/* Returns count of recognition warnings in warningsCountDataOut (if specified). */
IMAGO_CEXPORT int imagoRecognize(int *warningsCountDataOut);

/* Session-owned buffer, valid until the next Imago call on this session -
 * copy it out immediately (e.g. via QString::fromUtf8), same convention as
 * Indigo's own indigoMolfile()/indigoRxnfile(). No matching free function. */
IMAGO_CEXPORT const char *imagoGetMol();

#endif /* _imago_c_h */
