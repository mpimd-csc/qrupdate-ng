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

include Makeconf

help:
	@echo
	@echo "The following targets are available:"
	@echo "   make help    - displays this help"
	@echo "   make lib     - compiles a static library"
	@echo "   make solib   - compiles a dynamic library"
	@echo "   make test    - compiles and runs the testsuite"
	@echo "   make clean   - cleans up everything"
	@echo "   make install - installs everything"

lib:
	make -C src/ lib
solib:
	make -C src/ solib
test: lib
	make -C test/

clean:
	rm -f libqrupdate.a libqrupdate.so
	make -C src/ clean
	make -C test/ clean

install:
	make -C src/ install

install-shlib:
	make -C src/ install-shlib

install-staticlib:
	make -C src/ install-staticlib
