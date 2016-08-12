      ! adapted from MESA

      ! if vec contains decreasing values,
      ! returns 1 if val > vec(1); returns n if val <= vec(n)
      ! else returns k between 1 and n-1 such that
      !     vec(k) >= val > vec(k+1)

      ! if vec contains increasing values,
      ! returns 0 if val < vec(1); returns n if val >= vec(n)
      ! else returns k between 1 and n-1 such that
      !     vec(k) <= val < vec(k+1)

      integer function binary_search(n, vec, hint, val) result(loc)
         use const_def, only: dp
         integer, intent(in) :: n, hint
         real(dp), intent(in) :: val
         real(dp), intent(in) :: vec(:) ! (n)
         real(dp), parameter :: tiny = 1.0d-13
         integer :: first, last, mid
         
         if (n <= 1) then
            loc = n; return
         end if

         if (vec(n) < vec(1)) then ! decreasing values
         
            if (val > vec(1)) then
               loc = 0; return
            else if (abs(val - vec(1)) < tiny ) then
               loc = 1; return
            else if (val <= vec(n)) then
               loc = n; return
            end if
         
            if (hint > 0 .and. hint < n) then
               if (vec(hint) >= val .and. val > vec(hint+1)) then
                  loc = hint; return
               end if
            end if
            
            first = 1
            last = n-1
            loc = -1
            do while (first <= last)
               mid = (first + last)/2
               if (vec(mid) >= val) then
                  if (val > vec(mid+1)) then
                     loc = mid
                     exit
                  end if
                  first = mid + 1
               else
                  last = mid - 1
               end if
            end do

         else ! increasing values
         
            if (val < vec(1)) then
               loc = 0; return
            else if (abs(val - vec(1)) < tiny) then
               loc = 1; return
            else if (val >= vec(n)) then
               loc = n; return
            end if
         
            if (hint > 0 .and. hint < n) then
               if (vec(hint) <= val .and. val < vec(hint+1)) then
                  loc = hint; return
               end if
            end if
            
            first = 1
            last = n-1
            loc = -1
            do while (first <= last)
               mid = (first + last)/2
               if (vec(mid) <= val) then
                  if (val < vec(mid+1)) then
                     loc = mid
                     exit
                  end if
                  first = mid + 1
               else
                  last = mid - 1
               end if
            end do
               
         end if
         
      end function binary_search