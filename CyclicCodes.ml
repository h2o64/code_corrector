(* Polynomial type *)
type 'a poly = { deg : int ; coef : 'a array };;

(* Null polynomial *)
let null_poly () = {deg = -1 ; coef = [||]};;

(* Multiply with mod 2 *)
let mul_mod2 a b = (a * b) mod 2;;

(* Add with mod 2 *)
let add_mod2 a b = (a + b) mod 2;;

(* Sub with mod2 *)
let sub_mod2 a b = (a - b) mod 2;;

(* Array to polynomial *)
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

(* polynomials multiplication - Naive *)
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

(* Cyclotomic Poly - FIXME: URGENTLY *)
let cyclotomic_poly n = null_poly ();;

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

(* Substract polynomials *)
let rec sub_int_poly a b =
	(* Isolated cases *)
	if a.deg < 0 then fun_poly (~-) b
	else if b.deg < 0 then a
	else
		(* Actual substraction *)
		(if (a.deg > b.deg) then
			let ret = Array.copy a.coef in
			for i = 0 to b.deg do
				ret.(i)<-(sub_mod2 ret.(i) b.coef.(i));
			done;{deg = a.deg ; coef = ret}
		else if (b.deg > a.deg) then
			fun_poly (~-) (sub_int_poly b a)
		else
			(let ret = Array.copy a.coef in
			let ret_deg = ref (-1) in
			for i = 0 to a.deg do
				ret.(a.deg-i)<-(sub_mod2 ret.(a.deg-i) b.coef.(a.deg-i));
				if not (ret.(a.deg-i) = 0) && (!ret_deg < 0) then
					ret_deg := (a.deg-i);
			done;{deg = !ret_deg ; coef = ret}););;

(* Shift a polynomial *)
let shift_poly a n =
	let new_deg = a.deg + n in
	(* Copy to right index *)
	let new_coef = Array.make (a.deg+n+1) 0 in
	for i = 0 to a.deg do
		new_coef.(n+i)<-a.coef.(i);
	done;{deg = new_deg ; coef = new_coef};;

(* Print polynomial *)
let print_poly a =
	for i = 0 to a.deg do
		print_int a.coef.(a.deg-i);
		print_string "X^";
		print_int (a.deg-i);
		print_string " + ";
	done;(print_string "\n");;

(* Scalar multiplication for poly *)
let scalar_poly a l =
	(* Multiply a number by l *)
	let mul_aux x = mul_mod2 x l in
	if l = 0 then null_poly ()
	else { deg = a.deg ; coef = (Array.map mul_aux a.coef) };;

(* polynomial modulo - Naive *)
let rec mod_poly a n =
	if (a.deg < n.deg) then a
	else
		let shifted_poly = (shift_poly n (a.deg - n.deg)) in
		let adjusted_poly = scalar_poly shifted_poly (a.coef.(a.deg)) in
		mod_poly (sub_int_poly a adjusted_poly) n;;

(* Encode word *)
let encodeWord gen word_p n =
	(* Make the modulo polynomial *)
	let modulo = { deg = n ; coef = (Array.make (n+1) 0) } in
	modulo.coef.(0)<-(-1);
	modulo.coef.(n)<-1;
	(* Make a modulo *)
	let gen_mod = (mod_poly gen modulo) in
	(* Multiply the polynomial with generator *)
	let multiplied = (mul_poly word_p gen_mod add_mod2 mul_mod2 0) in
	multiplied;;
























