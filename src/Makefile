# Copyright (C) 2008, 2009  VZLU Prague, a.s., Czech Republic
#
# Author: Jaroslav Hajek <highegg@gmail.com>
#
# This file is part of qrupdate.
#
# qrupdate is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this software; see the file COPYING.  If not, see
# <http://www.gnu.org/licenses/>.
#
include ../Makeconf

SRC = caxcpy.f cch1dn.f cch1up.f cchdex.f cchinx.f cchshx.f cgqvec.f   \
cqhqr.f cqr1up.f cqrdec.f cqrder.f cqrinc.f cqrinr.f cqrot.f cqrqh.f   \
cqrshc.f cqrtv1.f dch1dn.f dch1up.f dchdex.f dchinx.f dchshx.f         \
dgqvec.f dqhqr.f dqr1up.f dqrdec.f dqrder.f dqrinc.f dqrinr.f dqrot.f  \
dqrqh.f dqrshc.f dqrtv1.f sch1dn.f sch1up.f schdex.f schinx.f schshx.f \
sgqvec.f sqhqr.f sqr1up.f sqrdec.f sqrder.f sqrinc.f sqrinr.f sqrot.f  \
sqrqh.f sqrshc.f sqrtv1.f zaxcpy.f zch1dn.f zch1up.f zchdex.f zchinx.f \
zchshx.f zgqvec.f zqhqr.f zqr1up.f zqrdec.f zqrder.f zqrinc.f zqrinr.f \
zqrot.f zqrqh.f zqrshc.f zqrtv1.f

OBJS = $(SRC:%.f=%.o)

lib: ../libqrupdate.a

solib: ../libqrupdate.so

../libqrupdate.a: $(OBJS)
	ar -cr $@ $(OBJS)

../libqrupdate.so: $(OBJS)
	$(FC) $(FFLAGS) -shared -o $@ -Wl,-soname=libqrupdate.so.$(MAJOR) $(OBJS) 

$(OBJS): %.o: %.f
	$(FC) $(FFLAGS) $(FPICFLAGS) -c $<

clean:
	rm -f $(OBJS)
