(* Use complex module *)
open Complex;;

(* Greatest Common Divisor *)
let rec gcd a b =
	match (a mod b) with
		| 0 -> b
		| r -> gcd b r;;

(* Polynome type *)
type 'a poly = { deg : int ; coef : 'a array };;

(* Null polynom *)
let null_poly () = {deg = -1 ; coef = [||]};;
 
(* Array to polynom *)
let poly_of_array arr zero =
	let getDeg a =
		let n = (Array.length a) in
		(* Find non-zero coeficient *)
		let ret = ref (-1) in
		for i = 0 to (n-1) do
			if a.(i) <> zero then ret := i
		done;!ret in
	{ deg = (getDeg arr);
		coef = arr };;

(* Polynoms multiplication - Naive *)
let mul_poly a b add_op mul_op zero =
	let ret_deg = if (a.deg < 0) || (b.deg < 0) then (-1) else (a.deg+b.deg) in
	let ret_coef = Array.make (if ret_deg > 0 then (ret_deg+1) else 0) zero in
	(* Make the coeficients *)
	if (ret_deg > 0) then
		(for i = 0 to a.deg do
			for j = 0 to b.deg do
				ret_coef.(i+j)<- add_op ret_coef.(i+j) (mul_op a.coef.(i) b.coef.(j));
			done;
		done);
	{deg = ret_deg ; coef = ret_coef};;

(* Compute cyclotomic polynom *)
let c_2 = { Complex.re = 2. ; Complex.im = 0.};;
let unitary_poly = poly_of_array [|Complex.one|] Complex.zero;;
let c_pi = { Complex.re = (4.0 *. atan 1.0) ; Complex.im = 0.};;
let cyclotomic_poly n =
	(* Quadruple complex multiplication and division *)
	let quad_mul a b c d e =
		Complex.div (Complex.mul (Complex.mul a b) (Complex.mul c d)) e in
	let ret = ref unitary_poly in
	let c_n = { Complex.re = (float_of_int n) ; Complex.im = 0. } in
	for k = 1 to n do
		if ((gcd k n) = 1) then
			let c_k = { Complex.re = (float_of_int k) ; Complex.im = 0. } in
			let cur = Complex.exp (quad_mul c_2 Complex.one Complex.i c_k c_n) in
			ret := mul_poly !ret (poly_of_array [|Complex.neg cur;Complex.one|] Complex.zero)
							Complex.add Complex.mul Complex.zero;
	done;
	(* Cleanup the polynom - FIXME: Slow and dirty *)
	let ret_new = { deg = !ret.deg ; coef = (Array.make (!ret.deg+1) 0)} in
	for i = 0 to !ret.deg do
		if (!ret.coef.(i).re > 0.5) then ret_new.coef.(i)<-1
		else if (!ret.coef.(i).re > 0.) then ret_new.coef.(i)<-0
		else if (!ret.coef.(i).re < (-0.5)) then ret_new.coef.(i)<-(-1)
		else ret_new.coef.(i)<-0;
	done;ret_new;;

(* Get all the possible codes *)
let getPolyCodes n =
	let rec getPolyCodes_aux a =
		if a = n then []
		else
			if not ((a mod n) = 0) then (cyclotomic_poly a)::(getPolyCodes_aux (a+1))
			else (getPolyCodes_aux (a+1)); in
	getPolyCodes_aux 1;;

(* Apply function on coefficients - Don't bother with degree *)
let fun_poly f a =
	{ deg = a.deg ; coef = (Array.map f a.coef) };;

(* Substract polynoms *)
let rec sub_int_poly a b =
	(* Isolated cases *)
	if a.deg < 0 then fun_poly (~-) b
	else if b.deg < 0 then a
	else
		(* Actual substraction *)
		(if (a.deg > b.deg) then
			let ret = Array.copy a.coef in
			for i = 0 to b.deg do
				ret.(i)<-(ret.(i) - b.coef.(i));
			done;{deg = a.deg ; coef = ret}
		else if (b.deg > a.deg) then
			fun_poly (~-) (sub_int_poly b a)
		else
			(let ret = Array.copy a.coef in
			let ret_deg = ref (-1) in
			for i = 0 to a.deg do
				ret.(a.deg-i)<-(ret.(a.deg-i) - b.coef.(a.deg-i));
				if not (ret.(a.deg-i) = 0) && (!ret_deg < 0) then
					ret_deg := (a.deg-i);
			done;{deg = !ret_deg ; coef = ret}););;

(* Shift a polynom *)
let shift_poly a n =
	let new_deg = a.deg + n in
	(* Copy to right index *)
	let new_coef = Array.make (a.deg+n+1) 0 in
	for i = 0 to a.deg do
		new_coef.(n+i)<-a.coef.(i);
	done;{deg = new_deg ; coef = new_coef};;

(* Print polynom *)
let print_poly a =
	for i = 0 to a.deg do
		print_int a.coef.(a.deg-i);
		print_string "X^";
		print_int (a.deg-i);
		print_string " + ";
	done;;

(* Scalar multiplication for poly *)
let scalar_poly a l =
	(* Multiply a number by l *)
	let mul_aux x = x * l in
	if l = 0 then null_poly ()
	else { deg = a.deg ; coef = (Array.map mul_aux a.coef) };;

(* Polynom modulo - Naive *)
let rec mod_poly a n =
	if (a.deg < n.deg) then a
	else
		let shifted_poly = (shift_poly n (a.deg - n.deg)) in
		let adjusted_poly = scalar_poly shifted_poly (a.coef.(a.deg)) in
		mod_poly (sub_int_poly a adjusted_poly) n;;

(* Encode word *)
let encodeWord gen word_p n =
	(* Make the modulo polynom *)
	let modulo = { deg = n ; coef = (Array.make (n+1) 0) } in
	modulo.coef.(0)<-(-1);
	modulo.coef.(n)<-1;
	(* Multiply the polynom with generator *)
	let multiplied = (mul_poly word_p gen (+) ( * ) 0) in
	(* Make a modulo *)
	mod_poly multiplied modulo;;



























