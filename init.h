void init_P(physics_grid *P,int problem);
void init_U( U_grid *U,physics_grid *P);
void init_F(F_grid *F,physics_grid *P);
physics_grid * create_physics_grid(void);
U_grid * create_U_grid(void);
F_grid * create_F_grid(void);
