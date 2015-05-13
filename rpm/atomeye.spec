Name:           atomeye
Version:        3
Release:        0
Url:            https://github.com/dmt4/atomeye
Summary:        AtomEye 3 atomistic configuration viewer
License:        Custom
Group:          Productivity/Scientific/Chemistry
Source:         %{name}-%{version}.tar.xz
BuildRoot:      %{_tmppath}/%{name}-%{version}-build
BuildRequires:  gcc-fortran
BuildRequires:  readline-devel
BuildRequires:  netcdf-devel
BuildRequires:  hdf5-devel
BuildRequires:  lapack-devel
BuildRequires:  gsl-devel
BuildRequires:  libX11-devel
BuildRequires:  libXext-devel
BuildRequires:  libXpm-devel
BuildRequires:  libpng16-compat-devel
BuildRequires:  libjpeg62-devel
BuildRequires:  fdupes


%description
The famous AtomEye viewer for atomistic simulations, large or small

This is a fork of James Kermode's fork of Ji Li's http://li.mit.edu/Archive/Graphics/A3/A3.html .


%prep
%setup -q


%build
%{__make} %{?_smp_mflags}


%install
install -D -s -m 755 bin/A %{buildroot}%{_bindir}/A
%fdupes -s $RPM_BUILD_ROOT


%clean
rm -rf %{buildroot}


%files
%defattr(-,root,root)
%doc README doc/
%{_bindir}/*


%changelog
