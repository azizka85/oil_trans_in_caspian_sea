__kernel void wind_induced_currents_davies85_variable_parameters_update_uniform_nu(
    float num,
    int ny, int nz,
    __global float* nu
) {
    int i = get_global_id(0);
    int j = get_global_id(1);
    int k = get_global_id(2);

    int p = j + i*ny;
    int id = k + p*nz;

    nu[id] = num;
}