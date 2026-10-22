! Export the written FAF records at full precision for independent regression.
program read_cap
  use mqc_gaussian
  implicit none
  type(mqc_gaussian_unformatted_matrix_file) :: f
  type(mqc_matrix) :: a
  character(:), allocatable :: path, record
  real(real64), allocatable :: row(:)
  integer :: i,j,u
  call mqc_get_command_argument(1,path)
  call mqc_get_command_argument(2,record)
  call f%getArray(record,a,fileName=path)
  open(newunit=u,file='record.txt',status='replace')
  write(u,*) size(a,1),size(a,2)
  allocate(row(size(a,2)))
  do i=1,size(a,1)
    do j=1,size(a,2)
      row(j)=a%at(i,j)
    end do
    write(u,'(*(ES25.16E3,1X))') row
  end do
  close(u)
  call Close_MatF(f%UnitNumber)
end program
