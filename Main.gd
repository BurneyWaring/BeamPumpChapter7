#From the Book:
#Create Compelling Science and Engineering Simulations Using the Godot Engine, Copyright 2024 Burney Waring, ThankGod Egbe, Lateef Kareem 
#Chapter 7

extends Node2D


var v1                       # point A
var v2                       # point B
var x3; var y3; var v3       # point C
var v4                       # point D

var vp0; var vp              # Initial and final vector representation of pitman
var vb0; var vb              # Initial and final vector representation of beam

var Lc;                      # Length of counterweight
var Lb2; var Lp2             # Length squared of beam and pitman respectively.
var smallRPM; var largeaRPM  # Revolution per minute for the small and large gear respectively.

var Rs        # Radius of the small gear
var Rl        # Radius of the large gear
var Rc        # radius of the gear in the counter weight.

var F        # matrix for functions in equation 7
var J        # Jacobian matrix
var delxy    # dx, dy

var res		#reference to equation function

func _ready():
	res = funcref(self, "_Equ")
	v1 = $counterweight.position
	v2 = $pitman.position
	v3 = $knob.position
	v4 = $beam.position
	
	vp0 = Vector2(v3-v2).normalized()
	vb0 = Vector2(v3-v4).normalized()
	
	smallRPM = 300
	
	Rs = 3
	Rl = 10 * Rs
	Rc = 20 * Rs
	
	Lc = sqrt((v2-v1).dot(v2-v1))
	Lb2 = (v3-v4).dot(v3-v4)
	Lp2 = (v3-v2).dot(v3-v2)

func _Equ(v):
	var vmv2 = v-v2
	var vmv4 = v-v4
	return Vector2(vmv2.dot(vmv2) - Lp2, vmv4.dot(vmv4) - Lb2)

func _process(delta):
	$display_motor_rotation.text = "Motor: "+ str(smallRPM)+" RPM"

	var anglechange = delta * smallRPM * 360.0/60.0
	$smallgear.rotation_degrees += anglechange
	
	anglechange = anglechange * Rs / Rl
	$largegear.rotation_degrees += anglechange
	
	anglechange = anglechange * Rl / Rc
	$counterweight.rotation_degrees += anglechange

	#SPM is stokes per minute and refers to the beam pump head motion, but
	# this is the same as the counterweight rotation per minute.
	#start with the angle change, degrees per frame.
	#anlgechange/delta is degrees per second.
	#anglechange/delta/360 is rotations per second
	#anglechange/delta/360*60 is rotations per minute, which is SPM	
	
	$display_SPM.text = "Pump: " + str(anglechange/delta/360*60)  + " SPM"
	
	if ($counterweight.rotation_degrees > 360.0):
		$counterweight.rotation_degrees -= 360.0
		
	var theta = $counterweight.rotation
	v2 = v1 + Lc*Vector2(cos(theta), sin(theta))
	$pitman.position = v2
	
	v3 = _FSolve(res, v3)
	
	$knob.position = v3
	vp = Vector2(v3-v2).normalized()
	vb = Vector2(v3-v4).normalized()

	$pitman.rotation = asin((vp0.cross(vp)))
	$beam.rotation = asin((vb0.cross(vb)))

func _on_motor_rpm_slider_value_changed(value):
	smallRPM = value
	
func _on_crankpin_radius_slider_value_changed(value):
	Lc = value # Replace with function body.
	
func _FSolve(fun, pt0):
	var pt = pt0
	var Fval = fun.call_func(pt)
	var error = Fval.length()
	while(error > 0.1):
		var Jac = _Jacobian(fun, pt)
		var dpt = _LinSolve(Fval, Jac)
		pt = pt-dpt
		Fval = fun.call_func(pt)
		error = Fval.length()
	return pt
	
#Solves Linear System
func _LinSolve(Fval, Jac):
	F = [Fval.x, Fval.y]
	J = [Jac[0].x, Jac[0].y, Jac[1].x, Jac[1].y]
	var det = J[0] * J[3] - J[1] * J[2]
	var dx = (J[3] * F[0] - J[2] * F[1]) / det
	var dy = (J[0] * F[1] - J[1] * F[0]) / det
	return Vector2(dx, dy)

#Computes the Jacobian of a function at point pt	
func _Jacobian(fun, pt):
	var f = fun.call_func(pt)
	var del = 1e-2
	var pdel
	var defx
	var defy
	pdel = pt
	pdel.x = pdel.x+del
	defx = (fun.call_func(pdel) - f)/del
	pdel = pt
	pdel.y = pdel.y+del
	defy = (fun.call_func(pdel) - f)/del
	return [defx, defy]
