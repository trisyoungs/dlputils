	! #####################
	! Periodic Table
	! #####################

	module ptable
	! Limits
        integer, parameter :: MAXEL = 300
	! Data
        character*8 :: elements(MAXEL)
	real*8 :: masses(MAXEL)
	integer :: nDefined = 0

	contains

	subroutine addElement(el, mass)
	implicit none
	character (len=*) :: el
	real*8 :: mass
	nDefined = nDefined + 1
	elements(nDefined) = el
	masses(nDefined) = mass
	end subroutine addElement

	subroutine initElements()
	implicit none

	! Populate initial element data?
	if (nDefined.eq.0) then
	  call addElement("H", 1.00798175d0)
	  call addElement("He", 4.002602d0)
	  call addElement("Li", 6.9675d0)
	  call addElement("Be", 9.0121831d0)
	  call addElement("B", 10.8135d0)
	  call addElement("C", 12.0106d0)
	  call addElement("N", 14.0069d0)
	  call addElement("O", 15.9994d0)
	  call addElement("F", 18.998403163d0)
	  call addElement("Ne", 20.1797d0)
	  call addElement("Na", 22.98976928d0)
	  call addElement("Mg", 24.3055d0)
	  call addElement("Al", 26.9815384d0)
	  call addElement("Si", 28.085d0)
	  call addElement("P", 30.973761998d0)
	  call addElement("S", 32.0675d0)
	  call addElement("Cl", 35.4515d0)
	  call addElement("Ar", 39.8775d0)
	  call addElement("K", 39.0983d0)
	  call addElement("Ca", 40.078d0)
	  call addElement("Sc", 44.955908d0)
	  call addElement("Ti", 47.867d0)
	  call addElement("V", 50.9415d0)
	  call addElement("Cr", 51.9961d0)
	  call addElement("Mn", 54.938043d0)
	  call addElement("Fe", 55.845d0)
	  call addElement("Co", 58.933194d0)
	  call addElement("Ni", 58.6934d0)
	  call addElement("Cu", 63.546d0)
	  call addElement("Zn", 65.38d0)
	  call addElement("Ga", 69.723d0)
	  call addElement("Ge", 72.630d0)
	  call addElement("As", 74.921595d0)
	  call addElement("Se", 78.971d0)
	  call addElement("Br", 79.904d0)
	  call addElement("Kr", 83.798d0)
	  call addElement("Rb", 85.4678d0)
	  call addElement("Sr", 87.62d0)
	  call addElement("Y", 88.90584d0)
	  call addElement("Zr", 91.224d0)
	  call addElement("Nb", 92.90637d0)
	  call addElement("Mo", 95.95d0)
	  call addElement("Tc", 98.0d0)
	  call addElement("Ru", 101.07d0)
	  call addElement("Rh", 102.90549d0)
	  call addElement("Pd", 106.42d0)
	  call addElement("Ag", 107.8682d0)
	  call addElement("Cd", 112.414d0)
	  call addElement("In", 114.818d0)
	  call addElement("Sn", 118.710d0)
	  call addElement("Sb", 121.760d0)
	  call addElement("Te", 127.60d0)
	  call addElement("I", 126.90447d0)
	  call addElement("Xe", 131.293d0)
	  call addElement("Cs", 132.90545196d0)
	  call addElement("Ba", 137.327d0)
	  call addElement("La", 138.90547d0)
	  call addElement("Ce", 140.116d0)
	  call addElement("Pr", 140.90766d0)
	  call addElement("Nd", 144.242d0)
	  call addElement("Pm", 145.0d0)
	  call addElement("Sm", 150.36d0)
	  call addElement("Eu", 151.964d0)
	  call addElement("Gd", 157.25d0)
	  call addElement("Tb", 158.925354d0)
	  call addElement("Dy", 162.500d0)
	  call addElement("Ho", 164.930328d0)
	  call addElement("Er", 167.259d0)
	  call addElement("Tm", 168.934218d0)
	  call addElement("Yb", 173.045d0)
	  call addElement("Lu", 174.9668d0)
	  call addElement("Hf", 178.49d0)
	  call addElement("Ta", 180.94788d0)
	  call addElement("W", 183.84d0)
	  call addElement("Re", 186.207d0)
	  call addElement("Os", 190.23d0)
	  call addElement("Ir", 192.217d0)
	  call addElement("Pt", 195.084d0)
	  call addElement("Au", 196.966570d0)
	  call addElement("Hg", 200.592d0)
	  call addElement("Tl", 204.383d0)
	  call addElement("Pb", 207.2d0)
	  call addElement("Bi", 208.98040d0)
	  call addElement("Po", 209.0d0)
	  call addElement("At", 210.0d0)
	  call addElement("Rn", 222.0d0)
	  call addElement("Fr", 223.0d0)
	  call addElement("Ra", 226.0d0)
	  call addElement("Ac", 227.0d0)
	  call addElement("Th", 232.0377d0)
	  call addElement("Pa", 231.03588d0)
	  call addElement("U", 238.02891d0)
	  call addElement("Np", 237.0d0)
	  call addElement("Pu", 244.0d0)
	  call addElement("Am", 243.0d0)
	  call addElement("Cm", 247.0d0)
	  call addElement("Bk", 247.0d0)
	  call addElement("Cf", 251.0d0)
	  call addElement("Es", 252.0d0)
	  call addElement("Fm", 257.0d0)
	  call addElement("Md", 258.0d0)
	  call addElement("No", 259.0d0)
	  call addElement("Lr", 262.0d0)
	  call addElement("Rf", 267.0d0)
	  call addElement("Db", 268.0d0)
	  call addElement("Sg", 269.0d0)
	  call addElement("Bh", 270.0d0)
	  call addElement("Hs", 269.0d0)
	  call addElement("Mt", 278.0d0)
	  call addElement("Ds", 281.0d0)
	  call addElement("Rg", 280.0d0)
	  call addElement("Cn", 285.0d0)
	  call addElement("Nh", 286.0d0)
	  call addElement("Fl", 289.0d0)
	  call addElement("Mc", 289.0d0)
	  call addElement("Lv", 293.0d0)
	  call addElement("Ts", 294.0d0)
	  call addElement("Og", 294.0d0)
	endif
	end subroutine initElements

	integer function findElement(el)
	implicit none
	character*8 :: el
	integer :: n

	! Ensure elements table is initialised
	call initElements()

	findElement = 0
	do n=1,nDefined
	  if (elements(n).eq.el) then
	    findElement = n
	    exit
	  end if
	end do
	end function findElement

	real*8 function getMass(el)
	implicit none
	character (len=*) :: el
	character*8 :: newEl
	integer :: i

	! Ensure elements table is initialised
	call initElements()

	! Try to find specified element in our table
	i = findElement(el)

	! If the element exists in the table already, return the mass
	if (i.ne.0) then
	  getMass = masses(i)
	  return
	end if

	! Not present, so request element name to assign to
	write(6,*) "Element ", el, "is unrecognised - please enter an element name from which to map its mass:"
	read(5,"(a4)") newEl
	i = findElement(newEl)
	if (i.ne.0) then
	  write(6,*) "Atom name ", el, "mapped to element ", newEl, "(mass = ", masses(i), ")"
	  ! Store the old unrecognised name in our map so if we find it again we know what mass it should be
	  call addElement(el, masses(i))
	  getMass = masses(i)
	  return
	end if

	! What to do?!?
	stop "Can't convert element to mass."

	end function getMass
	
	end module

