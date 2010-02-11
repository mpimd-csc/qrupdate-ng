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

LIBS= $(BLAS) $(LAPACK)
PROGS=tqr1up tqrinc tqrdec tqrshc tqrinr tqrder \
      tch1up tch1dn tchinx tchdex tchshx \
      tlu1up tlup1up

OUTS=$(PROGS:%=%.out)

tests: $(OUTS)
	./report_results $(OUTS)

$(OUTS): %.out: %
	echo > $@
	./$< | tee $@

$(PROGS): % : %.f utils.o ../libqrupdate.a
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

utils.o: utils.f
	$(FC) $(FFLAGS) -c $<

../libqrupdate.a:
	make -C ../ lib
	
clean: 
	rm -f *.o $(PROGS) $(OUTS)
