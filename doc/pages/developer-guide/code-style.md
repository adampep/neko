# Code style {#code-style}

\tableofcontents

Fortran code should conform to the Fortran 2008 standard and should use an
indentation level of 2, except for the extra indentation within `do` `if`,
`select` or `where` statements and for each level inside a structure e.g.
`type`, `interface`, where the indentation level is 3, continuation statements,
which should be indented by 5 and a 0 indentation is
used for `module` or `contains` (except for `contains` inside a derived type,
where a single indentation level is used). It should be noted that no
indentation should be applied to the `case` statements contained in a 
`select case` block.

These are the default rules in Emacs' Fortran mode, an example is given below,
additional information on the Emacs' Fortran mode can be found at
[https://emacsdocs.org](https://emacsdocs.org/docs/emacs/Fortran) and the full
formatting script for emacs can be found at the [Emacs' GitHub](https://github.com/emacs-mirror/emacs/blob/master/lisp/progmodes/f90.el).

~~~~~~~~~~~~~~~{.f90}
module example
  use mod
  implicit none

  type :: derived_t
     integer :: x
   contains
     procedure, pass(this) :: bar
  end type derived_t

contains

  subroutine foo(x, y)
    integer, intent(in) :: x
    integer, intent(inout) :: y
    real(kind=rp) :: value

    do i = 1, 10
       ...
    end do

    if (x .lt. y) then
       ... This statement is very long &
            continues on the next line
    end if
    
    select case (x)
    case (1)
       ...
    case (2)
       ...
    end select

  end subroutine foo

end module example
~~~~~~~~~~~~~~~

## Data types
For portability reasons, it is essential to only use data type kinds defined in
num_types.f90 and avoid legacy constructs like `real*8` or `integer(kind=8)`

Floating-point numbers should be declared using `real` with the kind `rp`, which
is the configured working precision for Neko (defaults to double). If single,
double or quad precision is explicitly needed, use the kinds `sp`, `dp` or `qp`,
respectively. For 16, 32, or 64-bit integers, Neko has defined the kinds ` i2`,
`i4` or `i8`, respectively; however, for standard integers, it is perfectly fine
to use a plain `integer`.

## Enforcing rules

Submitting a pull request to Neko requires that the code adheres to the linting
rules and to the formatting. The linting will require that the code does not
reduce the linting score of the repository, while the formatting should be
completed for any modified files.

In addition to the indentation specified above, the following rules are
enforced:

- Lines may not exceed 80 characters.
- Do loop specification must have spaces `do i = 1, 10`.
- Logical operators must have spaces around them `a .eq. b`.
- The separator `::` must have a spaces around it when declaring variables.
- Punctuations must have spaces after them `foo(b, c)`.
- Context blocks must have a space before the parenthesis `if (a .eq. b)`.
- Usage of OpenMP should be prepended with `!$`.
- Indentation should be done with spaces, not tabs.
- Use new syntax `type(kind)` instead of `type*8`.
- Comment operator `!` must have a space before and after them `! This is a comment`.
- Lines may not be terminated with a semicolon.
- End statement should have a context block `end if`.
- End statement should have a space before context specification `end if`.
- Assignment operators must have spaces around them `a = b`.
- Trailing white spaces is not allowed.
- Double spaces are not allowed.
- Precision of real numbers should be specified using `sp`, `dp` or `qp`.
- Array declaration should use brackets instead of parentheses.
- Should use `use mpi_f08` instead (or `use mpi` if not available).
- Bare stop statement not allowed.
- The use of bare `exit` statement is not allowed.
- The use of the `goto` statement is not allowed.
- The use of the `pause` statement is not allowed.
- The use of the `include` statement is not allowed.

However, there are some exceptions to these rules:

- The separator `::` may not have spaces around it in the case of a type
  declaration.
- Spaces after the comma in a list of arguments or array indices may be omitted
  for single letter variable names and 2 digit numbers. For example, `foo(a,b)`
  and `foo(a,10)` are allowed, but `foo(a,bb)` and `foo(a,100)` are not.
- Spaces after the comma are not required in format specifiers.
- Spaces around the `=` operator are not required in type declarations.

### Tools

In order to simplify compliance to the formatting and linting rules the
following tools can be used:

- [flint](https://github.com/marshallward/flint) to test the linting compliance.
- [findent](https://github.com/wvermin/findent) to apply the formatting rules.

We have scripts available under the `contrib/lint_format` directory that can be
used to lint and format the code. The scripts are `lint.sh` and `format.sh`
respectively.

#### flint

The linting rules are enforced by the
[flint](https://github.com/marshallward/flint) tool. The rules are defined in the
`flinter_rc.yml` file in the root of the repository. 

Please note, newer versions of flint fails to execute for some of our large
files (htable.f90 and stack.f90).

One way to install flint is through pip:
```sh
pip install nobvisual==0.2.0 flinter==0.4.0
```

To test your code against the linting rules, you can run the following command:

```sh
flint score -r flinter_rc.yml <file>.f90
```
The whole src directory can be checked with:
```sh
flint score -r flinter_rc.yml src/
```

#### findent

This tool can be used to enforce these rules by assigning the following options.
The documentation of `findent` provide details for emacs, vim and gedit. For
VSCode, the [Modern Fortran](https://marketplace.visualstudio.com/items?itemName=fortran-lang.linter-gfortran)
extension provides an integration.

```sh
pip install findent
findent -Rr -i2 -d3 -f3 -s3 -w3 -t3 -j3 -k- --ws_remred --indent_ampersand --openmp=0 < file.f90 > file.f90.tmp
mv file.f90.tmp file.f90
```
