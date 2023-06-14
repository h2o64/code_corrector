(* Greatest Common Divisor *)
let rec gcd a b =
	match (a mod b) with
		| 0 -> b
		| r -> gcd b r;;

(* Prime Factorization - Pollard's Rho Algorithm *)
let pollard_factor a =
	(* Hardcoded initial parameters *)
	let x_fixed = ref 2 in
	let cycle_size = ref 2 in
	let x = ref 2 in
	let factor = ref 1 in
	while (!factor = 1) do
		let count = ref 1 in
		while (!count <= !cycle_size) && (!factor <= 1) do	
			(* Use P(X) = x^2 + 1 *)
			x := ((!x * !x) + 1) mod a;
			factor := gcd (abs (!x - !x_fixed)) a;
			count := !count + 1;
		done;
		cycle_size := !cycle_size * 2;
		x_fixed := !x;
	done;!factor;;

(* Actual prime decomposition *)
let rec decomposition factorization n =
	if (n = 1) || (n < 0) then []
	else
		(* Do the decomposition *)
		(let factor = (factorization n) in
		factor::(decomposition factorization (n /factor)));;

(* Count multiplicity in list - lists are small *)
let rec multiplicity l =
	(* Let's max it out to 20 *)
	let seen = Array.make 20 0 in
	let seen_cur = ref 0 in
	let rec multiplicity_aux f (x,m) = match f with
		| [] -> (x,m)
		| h::t -> if (h = x) then multiplicity_aux t (x,(m+1))
							else multiplicity_aux t (x,m) in
	match l with
		| [] -> []
		| h::t ->
			if not (Array.mem h seen) then
				(seen.(!seen_cur)<-h;
				seen_cur := !seen_cur + 1;
				(multiplicity_aux l (h,0))::(multiplicity t))
			else (multiplicity t);;

(* Is there a item with multiplicity higher than one *)
let rec is_higher_than_one_mul l = match l with
	| [] -> false
	| (_,a)::t -> if a > 1 then true else (is_higher_than_one_mul t);;

(* Multiply a whole list *)
let rec mul_list l = match l with
	| [] -> 1
	| h::t -> h*(mul_list t);;

(* Let first element of a tuple list - apply a function to each *)
let rec mul_tuple_list l f = match l with
	| [] -> 1
	| (a,_)::t -> ((f a)*(mul_tuple_list t f));;

(* Batman method to evaluate coefficient sizes - UNFINISHED *)
let rec batman_limit j k l =
	if k = j then 1
	else
		match l with
			| [] -> failwith "batman_limit: Fail";
			| h::t -> BitsTools.pow h ((BitsTools.pow 2 (j-k-2))-1) * (batman_limit j (k+1) t);;

(* Make list superset *)
let rec superset = function
  | [] -> [[]]
  | x :: xs -> 
     let ps = superset xs in
     ps @ List.map (fun ss -> x :: ss) ps;;

(* Get the coefficients - Based on work from Andrew Arnold and Michael Monagan *)
(* Based of Sage's implementation *)
let rec cyclotomic_coef n =
	(* Get prime factorization of n *)
	let factors = (decomposition pollard_factor n) in
	let factors_m = multiplicity factors in
  (* If there are primes that occur in the factorization with multiplicity
     greater than one we use the fact that Phi_ar(x) = Phi_r(x^a) when all
     primes dividing a divide r. *)
	if (is_higher_than_one_mul factors_m) then
		(print_string "HER";
		let f x = x in
		let g x = (x-1) in
		let rad = mul_tuple_list factors_m f in
		let rad_coeffs = cyclotomic_coef rad in
		let pw = n / rad in
		let l = Array.make 0 (1+pw*(mul_tuple_list factors_m g)) in
		for i = 0 to ((Array.length rad_coeffs)-1) do
			l.(i*pw)<-rad_coeffs.(i);
		done;
		rad_coeffs)
	(* If n is prime, it's easy *)
	else if ((List.length factors_m) = 1) then
		(Array.make n 1)
	else
		begin
		(* The following bounds are from Michael Monagan:
       For all n < 169,828,113, the height of Phi_n(x) is less than 60 bits.
       At n = 169828113, we get a height of 31484567640915734951 which is 65 bits,
       so we fail
       For n=10163195, the height of Phi_n(x) is 1376877780831,  40.32 bits.
       For n<10163195, the height of Phi_n(x) is <= 74989473, 26.16 bits. *)
		let fits_long_limit = 10163195 in
		(* Evaluate the coefficient sizes with bateman method *)
		(* let batman_bound = batman_limit ((List.length factors)-2) 0 factors in
		if batman_bound > max_int then failwith "cyclotomic_coef: Coefficients too big"; *)
		if n > fits_long_limit then failwith "cyclotomic_coef: Coefficients too big";
		(* Sort the factors everytime *)
		let prime_subsets = superset (List.stable_sort compare factors) in
		(* Get the maximum degree *)
		let rec max_deg l d = match l with
			| [] -> d
			| s::t ->
				if ((List.length s) mod 2) = 0 then
					let tmp_prod = mul_list s in
					max_deg t (d + (n / tmp_prod))
				else max_deg t d in
		let max_deg_p = max_deg prime_subsets 0 in
		(* Use palyndrum properties of cyclotomic polynomials *)
		let coeffs = Array.make (max_deg_p+2) 0 in
		(* Set constant coef to 1 *)
		coeffs.(0)<-1;
		let offset = ref 0 in
		let deg = ref 0 in
		print_string "HERE1";
		let rec first_side l = match l with
			| [] -> ()
			| s::t ->
				if (((List.length s) mod 2) = 0) then
					(let d = mul_list s in
					let dd = n / d in
					for k = dd to (!deg+dd) do
						coeffs.(k) <- coeffs.(k) - coeffs.(k-dd);
					done;
					deg := !deg + dd;);
				first_side t in
		first_side prime_subsets;
		print_string "HERE2";
		let rec second_side l = match l with
			| [] -> ()
			| s::t ->
				if (((List.length s) mod 2) = 1) then
					(let d = mul_list s in
					let dd = n / d in
					for k = !deg-dd+1 to !deg do
						coeffs.(k) <- (-1) * coeffs.(k);
					done;
					for k = !offset to !deg-dd do
						coeffs.(k) <- coeffs.(k+dd) - coeffs.(k);
					done;
					offset := !offset + dd;);
				second_side t in
		second_side (List.rev prime_subsets);
		print_string "HERE3\n";
		print_int !deg;
		print_string "\n";
		print_int max_deg_p;
		print_string "\n";
		print_int !offset;
		let ret = Array.make (!deg+1) 0 in
		let k = ref (!offset-1) in
		while (!k < !deg+1) do
			ret.(!k - !offset + 1)<-coeffs.(!k);
			k := !k + 1;
		done;
		ret;
		end;;
















