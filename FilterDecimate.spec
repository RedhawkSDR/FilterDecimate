#
#This file is protected by Copyright. Please refer to the COPYRIGHT file
# distributed with this source distribution.

# This file is part of FilterDecimate.
#
# FilterDecimate is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.
#
# FilterDecimate is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
# for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see http://www.gnu.org/licenses/.
#

# By default, the RPM will install to the standard REDHAWK SDR root location (/var/redhawk/sdr)
# You can override this at install time using --prefix /new/sdr/root when invoking rpm (preferred method, if you must)
%{!?_sdrroot: %define _sdrroot /var/redhawk/sdr}
%define _prefix %{_sdrroot}
Prefix: %{_prefix}

# Point install paths to locations within our target SDR root
%define _sysconfdir    %{_prefix}/etc
%define _localstatedir %{_prefix}/var
%define _mandir        %{_prefix}/man
%define _infodir       %{_prefix}/info

Name: FilterDecimate
Summary: Component %{name}
Version: 1.0.0
Release: 1
License: None
Group: REDHAWK/Components
Source: %{name}-%{version}.tar.gz
BuildRoot: %{_tmppath}/%{name}-root

Requires: redhawk >= 1.9
BuildRequires: redhawk-devel >= 1.9
BuildRequires: autoconf automake libtool

# Interface requirements
Requires: bulkioInterfaces
BuildRequires: bulkioInterfaces

%description
Component %{name}


%prep
%setup


%build
# Implementation cpp
pushd cpp
./reconf
%define _bindir %{_prefix}/dom/components/FilterDecimate/cpp
%configure
make %{?_smp_mflags}
popd
# Implementation cpp_arm
pushd cpp_arm
./reconf
%define _bindir %{_prefix}/dom/components/FilterDecimate/cpp_arm
%configure
make %{?_smp_mflags}
popd


%install
rm -rf $RPM_BUILD_ROOT
# Implementation cpp
pushd cpp
%define _bindir %{_prefix}/dom/components/FilterDecimate/cpp
make install DESTDIR=$RPM_BUILD_ROOT
popd
# Implementation cpp_arm
pushd cpp_arm
%define _bindir %{_prefix}/dom/components/FilterDecimate/cpp_arm
make install DESTDIR=$RPM_BUILD_ROOT
popd


%clean
rm -rf $RPM_BUILD_ROOT


%files
%defattr(-,redhawk,redhawk)
%dir %{_prefix}/dom/components/%{name}
%{_prefix}/dom/components/%{name}/FilterDecimate.scd.xml
%{_prefix}/dom/components/%{name}/FilterDecimate.prf.xml
%{_prefix}/dom/components/%{name}/FilterDecimate.spd.xml
%{_prefix}/dom/components/%{name}/cpp
%{_prefix}/dom/components/%{name}/cpp_arm

