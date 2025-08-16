__kernel void wind_induced_currents_davies85_variable_parameters_update_uniform_q(
    float qm,
    int ny,
    __global float* q
) {
    int i = get_global_id(0);
    int j = get_global_id(1);

    int k = j + i*ny;

    q[k] = qm;
}