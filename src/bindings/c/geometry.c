
#include "geometry.h"

#include "string.h"

extern void add_matching_mesh_c(const char *umbilical_elem_option, int *umbilical_elem_option_len, int umbilical_element_numbers[], int *umbilical_element_numbers_len);
void append_units_c();
void create_anastomosis_c(int *node_in, int *node_out, double *anast_radius);
void define_1d_element_geometry_c(const char *ELEMFILE, int *filename_len);
void define_1d_element_placenta_c(const char *ELEMFILE, int *filename_len, int *anastomosis_elem_in);
void define_node_geometry_c(const char *NODEFILE, int *filename_len);
void define_rad_from_file_c(const char *FIELDFILE, int *filename_len, const char *order_system, 
                            int *order_system_len, double *s_ratio);
void define_rad_from_geom_c(const char *order_system, int *order_system_len, double *control_param,
                            const char *start_from, int *start_from_len, double *start_rad,
                            const char *group_type, int *group_type_len, const char *group_options, int *group_options_len);
void define_ven_rad_from_art_c(const char *FILENAME, int *filename_len, double *factor);
void define_capillary_model_c(int *define_convolutes, int *define_generations, int *define_parallel,const char *DEFINE_MODEL,int *define_model_len);
void element_connectivity_1d_c();
void evaluate_ordering_c();
void update_radius_by_order_c(int *order,double *new_radius,const char *UPDATE_TYPE, int *update_type_len, double *control_parameter);
void update_1d_elem_field_c(int *ne_field, int *elem_number, double *value);


void add_matching_mesh(const char *umbilical_elem_option, int umbilical_element_numbers_len, int umbilical_element_numbers[])
{
  int umbilical_elem_option_len = strlen(umbilical_elem_option);
  add_matching_mesh_c(umbilical_elem_option, &umbilical_elem_option_len, umbilical_element_numbers, &umbilical_element_numbers_len);
}

void append_units()
{
  append_units_c();
}


void create_anastomosis(int node_in, int node_out, double anast_radius)
{
  create_anastomosis_c(&node_in, &node_out, &anast_radius);
}

void define_1d_element_geometry(const char *ELEMFILE)
{
  int filename_len = strlen(ELEMFILE);
  define_1d_element_geometry_c(ELEMFILE, &filename_len);
}

void define_1d_element_placenta(const char *ELEMFILE,int anastomosis_elem_in)
{
  int filename_len = strlen(ELEMFILE);
  define_1d_element_placenta_c(ELEMFILE, &filename_len, &anastomosis_elem_in);
}


void define_node_geometry(const char *NODEFILE)
{
  int filename_len = strlen(NODEFILE);
  define_node_geometry_c(NODEFILE, &filename_len);
}

void define_rad_from_file(const char *FIELDFILE,const char *order_system,double s_ratio)
{
  int filename_len = strlen(FIELDFILE);
  int order_system_len = strlen(order_system);
  define_rad_from_file_c(FIELDFILE, &filename_len, order_system, &order_system_len, &s_ratio);
}

void define_rad_from_geom(const char *order_system, double control_param, const char *start_from,
                          double start_rad, const char*group_type, const char *group_options)
{
  int order_system_len = strlen(order_system);
  int start_from_len = strlen(start_from);
  int group_type_len = strlen(group_type);
  int group_options_len = strlen(group_options);
  define_rad_from_geom_c(order_system, &order_system_len, &control_param, start_from, &start_from_len, &start_rad,
                         group_type, &group_type_len, group_options, &group_options_len);

}

void define_ven_rad_from_art(const char *FILENAME,double factor)
{
  int filename_len = strlen(FILENAME);
  define_ven_rad_from_art_c(FILENAME, &filename_len, &factor);
}

void define_capillary_model(int define_convolutes, int define_generations, int define_parallel,const char *DEFINE_MODEL)
{
  int define_model_len = strlen(DEFINE_MODEL);
  define_capillary_model_c(&define_convolutes, &define_generations, &define_parallel, DEFINE_MODEL, &define_model_len);
}

void element_connectivity_1d()
{
  element_connectivity_1d_c();
}

void update_radius_by_order(int order,double new_radius,const char *UPDATE_TYPE, double control_parameter)
{
  int update_type_len = strlen(UPDATE_TYPE);
  update_radius_by_order_c(&order, &new_radius, UPDATE_TYPE, &update_type_len, &control_parameter);
}

void evaluate_ordering()
{
  evaluate_ordering_c();
}

void update_1d_elem_field(int ne_field, int elem_number, double value)
{
  update_1d_elem_field_c(&ne_field, &elem_number, &value);
}
