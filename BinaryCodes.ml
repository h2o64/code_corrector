module BinaryCodes :
	sig
		val even_code : int -> int -> int
		val even_code_improved : int -> int
		val repeat_code : int -> int -> int
	end =
	struct

	(* Even code *)
	(* If the weight is even, add a 0 else add 1 *)
	let even_code a size =
		let w = BitsTools.hamming_weight a size in
		let shiffted_a = (a lsl 1) in 
		if (w mod 2) = 0 then shiffted_a
		else (shiffted_a + 1);;

	(* Even code - v2 *)
	(* The number is split in may 2-long numbers *)
	let even_code_improved a size =
		let rec even_code_improved_aux cur ret iter =
			if cur = 0 then ret
			else
				even_code_improved_aux (cur lsr 2)
			                         (ret + ((even_code (cur land 3) size) lsl (iter*3)))
			                         (iter+1) in
		even_code_improved_aux a 0 0;;

	(* Repeat code *)
	(* Each bit gets repeated n times *)
	let repeat_code a n =
		let i = ref 0 in
		let cur = ref a in
		let ret = ref 0 in
		let pow_n = ((BitsTools.pow 2 n)-1) in (* FIXME: Make a hardcoded table *)
		while not (!cur = 0) do
			if ((!cur land 1) = 1) then
				ret := (!ret + (pow_n lsl ((!i)*n)));
			cur := (!cur lsr 1);
			i := !i + 1
		done;!ret;;

	end
