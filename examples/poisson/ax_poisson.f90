module ax_poisson
  use ax_product
  use mxm_wrapper
  use mesh
  use space
  use coefs
  use math
  use num_types
  use utils
  implicit none

  type, public, extends(ax_t) :: ax_poisson_t
  contains
    procedure, nopass :: compute => ax_poisson_compute
    procedure, pass(this) :: compute_vector => ax_poisson_compute_vector
  end type ax_poisson_t

contains 

  subroutine ax_poisson_compute(w, u, coef, msh, Xh)
    type(mesh_t), intent(inout) :: msh
    type(space_t), intent(inout) :: Xh
    type(coef_t), intent(inout) :: coef
    real(kind=rp), intent(inout) :: w(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(inout) :: u(Xh%lx, Xh%ly, Xh%lz, msh%nelv)

    real(kind=rp) :: ur(Xh%lx**3)
    real(kind=rp) :: us(Xh%lx**3)
    real(kind=rp) :: ut(Xh%lx**3)
    real(kind=rp) :: wk(Xh%lx**3)
    integer :: e

    do e = 1, msh%nelv
       call ax_e(w(1,1,1,e), u(1,1,1,e), &
            coef%g11(1,1,1,e), coef%g22(1,1,1,e), coef%g33(1,1,1,e), &
            ur, us, ut, wk, Xh%lx, Xh%dx, Xh%dxt)
    end do

  end subroutine ax_poisson_compute
  
  subroutine ax_e(w, u, g11, g22, g33, ur, us, ut, wk, lx, D, Dt)
    integer, intent(inout) :: lx
    real(kind=rp), intent(inout) :: w(lx**3)
    real(kind=rp), intent(inout) :: u(lx**3)
    real(kind=rp), intent(inout) :: g11(lx**3)
    real(kind=rp), intent(inout) :: g22(lx**3)
    real(kind=rp), intent(inout) :: g33(lx**3)
    real(kind=rp), intent(inout) :: ur(lx**3)
    real(kind=rp), intent(inout) :: us(lx**3)
    real(kind=rp), intent(inout) :: ut(lx**3)
    real(kind=rp), intent(inout) :: wk(lx**3)
    real(kind=rp), intent(inout) :: D(lx, lx)
    real(kind=rp), intent(inout) :: Dt(lx, lx)
    real(kind=rp) :: wr, ws, wt
    integer :: i, n
  
    n = lx - 1
    call local_grad3(ur, us, ut, u, n, D, Dt)

    do i = 1, lx**3
       wr = g11(i)*ur(i) 
       ws = g22(i)*us(i) 
       wt = g33(i)*ut(i) 
       ur(i) = wr
       us(i) = ws
       ut(i) = wt
    end do

    call local_grad3_t(w, ur, us, ut, n, D, Dt, wk)
  
  end subroutine ax_e
  
  subroutine local_grad3(ur, us, ut, u, n, D, Dt)
    integer, intent(inout) :: n    
    real(kind=rp), intent(inout) :: ur(0:n, 0:n, 0:n)
    real(kind=rp), intent(inout) :: us(0:n, 0:n, 0:n)
    real(kind=rp), intent(inout) :: ut(0:n, 0:n, 0:n)
    real(kind=rp), intent(inout) :: u(0:n, 0:n, 0:n)
    real(kind=rp), intent(inout) :: D(0:n, 0:n)
    real(kind=rp), intent(inout) :: Dt(0:n, 0:n)
    integer :: m1, m2, k
  
    m1 = n + 1
    m2 = m1*m1
  
    call mxm(D ,m1,u,m1,ur,m2)
    do k = 0, n
       call mxm(u(0,0,k),m1,Dt,m1,us(0,0,k),m1)
    end do
    call mxm(u,m2,Dt,m1,ut,m1)
    
  end subroutine local_grad3
  
  subroutine local_grad3_t(u,ur,us,ut,N,D,Dt,w)
    integer, intent(inout) :: n  
    real(kind=rp), intent(inout) :: ur(0:n, 0:n, 0:n)
    real(kind=rp), intent(inout) :: us(0:n, 0:n, 0:n)
    real(kind=rp), intent(inout) :: ut(0:n, 0:n, 0:n)
    real(kind=rp), intent(inout) :: u(0:n, 0:n, 0:n)
    real(kind=rp), intent(inout) :: D(0:n, 0:n)
    real(kind=rp), intent(inout) :: Dt(0:n, 0:n)
    real(kind=rp), intent(inout) :: w(0:n, 0:n, 0:n)
    integer :: m1, m2, m3, k
  
    m1 = n + 1
    m2 = m1**2
    m3 = m1**3
    
    call mxm(Dt,m1,ur,m1,u,m2)
    
    do k = 0, N
       call mxm(us(0,0,k),m1,D ,m1,w(0,0,k),m1)
    end do
    call add2(u,w,m3)
    
    call mxm(ut,m2,D ,m1,w,m1)
    call add2(u,w,m3)
    
  end subroutine local_grad3_t

  subroutine ax_poisson_compute_vector(this, au, av, aw, u, v, w, coef, msh, Xh)
    class(ax_poisson_t), intent(in) :: this
    type(space_t), intent(inout) :: Xh
    type(mesh_t), intent(inout) :: msh
    type(coef_t), intent(inout) :: coef
    real(kind=rp), intent(inout) :: au(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(inout) :: av(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(inout) :: aw(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(inout) :: u(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(inout) :: v(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(inout) :: w(Xh%lx, Xh%ly, Xh%lz, msh%nelv)

    call neko_error('Not in use')

  end subroutine ax_poisson_compute_vector

end module ax_poisson
