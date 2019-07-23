//***************************************************************
//***************************************************************

void ODE::setup(Malaria_state &y0, Parms& parms0, const double t0, const double h0) {
	step_type = from_map_int("step_type", 1, 3, 1);
	t = t0;
	step_size = h0;
	y = &y0;
	parms = &parms0;
	resize();
}
void ODE::resize() {
	ytemp = *y;
	k1 = *y;
	k2 = *y;
	k3 = *y;
	k4 = *y;
}
void ODE::step(const double total_step){
	const double end_time = t + total_step;
	const int m = ceil(total_step/step_size);
	const double h = total_step/m;
	for(int j=0; j<m; j++){
		parms->update_history(*y, t);
		sub_step(h);
	}
	t = end_time;
	parms->update_history(*y, t);
}

void ODE::sub_step(const double h){
	const int n = parms->run_type==2 ? 6*num_mosq : y->size_now();
	const int s = parms->run_type==2 ? 3 : step_type;
	if (s==1) {
		parms->calc_dydt(t, *y, k1);
		for (int i = 0; i<n; i++)
			(*y)[i] = max((*y)[i]+h*k1[i], 0.0);
	}
	else if (s==2) {
		parms->calc_dydt(t, *y, k1);
		for (int i = 0; i<n; i++)
			ytemp[i] = (*y)[i]+h*k1[i];
		parms->calc_dydt(t+h, ytemp, k2);
		for (int i = 0; i<n; i++)
			(*y)[i] = max((*y)[i]+h*0.5*(k1[i] + k2[i]), 0.0);
	}
	else if (s==3) {
		parms->calc_dydt(t, *y, k1);
		for (int i = 0; i<n; i++)
			ytemp[i] = (*y)[i]+h*0.5*k1[i];
		parms->calc_dydt(t+h*0.5, ytemp, k2);
		for (int i = 0; i<n; i++)
			ytemp[i] = (*y)[i]+h*0.5*k2[i];
		parms->calc_dydt(t+h*0.5, ytemp, k3);
		for (int i = 0; i<n; i++)
			ytemp[i] = (*y)[i]+h*k3[i];
		parms->calc_dydt(t+h, ytemp, k4);
		for (int i = 0; i<n; i++)
			(*y)[i] = max((*y)[i]+h/6.0*(k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]), 0.0);
	}
	t += h;
}
//***************************************************************
//***************************************************************