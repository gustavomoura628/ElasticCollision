#include <SDL2/SDL.h> 
//#include <SDL2/SDL_image.h> 
#include <SDL2/SDL_timer.h> 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <math.h>
const int SCREEN_WIDTH = 1366;
const int SCREEN_HEIGHT = 768;
const int FPS = 60;
const float deltaTime = 1;
const float PI = 3.14159265358979;
void CSDL_DrawCircle(SDL_Renderer* renderer, int x, int y, int r){
	int i,j;
	/*
	for(i=0;i<r;i++){
		for(j=-sqrt(r*r - i*i);j<+sqrt(r*r - i*i);j++){
			SDL_RenderDrawPoint(renderer,x+j,y+i);
			SDL_RenderDrawPoint(renderer,x+j,y-i);
		}
	}
	*/
	j=0;
	i=r;
	int k;
	int f=0;
	for(j=0;i>=j;j++){
		f=j*j+i*i-r*r;
		if(f>0)i--;
		SDL_RenderDrawLine(renderer, x-j, y+i, x+j, y+i);
		SDL_RenderDrawLine(renderer, x-j, y-i, x+j, y-i);
		SDL_RenderDrawLine(renderer, x-i, y+j, x+i, y+j);
		SDL_RenderDrawLine(renderer, x-i, y-j, x+i, y-j);
		/*
		for(k=-j;k<=j;k++){
			SDL_RenderDrawPoint(renderer,x+k,y+i);
			SDL_RenderDrawPoint(renderer,x+k,y-i);
		}
		for(k=-i;k<=i;k++){
			SDL_RenderDrawPoint(renderer, x+k,y+j);
			SDL_RenderDrawPoint(renderer, x+k,y-j);
		}
		*/
	}
}


typedef struct vector{
	float x;
	float y;
}Vector;

typedef struct particle{
	float mass;
	float radius;
	Vector pos;
	Vector vel;
	int color_r;
	int color_g;
	int color_b;
}Particle;

void vector_sum(Vector * v0, Vector * v1, Vector * vr){
	vr->x = v0->x + v1->x;
	vr->y = v0->y + v1->y;
}
void vector_subtract(Vector * v0, Vector * v1, Vector * vr){
	vr->x = v0->x - v1->x;
	vr->y = v0->y - v1->y;
}
float vector_module(Vector * v){
	return sqrt(v->x*v->x + v->y*v->y);
}
void vector_projection(Vector * a, Vector * b, Vector * r){
	float s = (a->x * b->x + a->y * b->y)/(b->x*b->x+b->y*b->y);
	r->x = s*b->x;
	r->y = s*b->y;
}
void vector_rotate(Vector * v, Vector * r, float angle){
	r->x = v->x * cos(angle) - v->y*sin(angle);
	r->y = v->x * sin(angle) + v->y*cos(angle);
}
void vector_multiplyReal(Vector * v, Vector * r, float value){
	r->x=value*v->x;
	r->y=value*v->y;
}
float vector_dotProduct(Vector * v0, Vector * v1){
	return (v0->x * v1->x)+(v0->y * v1->y);
}
void vector_normalize(Vector * v, Vector * r){
	r->x = v->x/vector_module(v);
	r->y = v->y/vector_module(v);
}
void particle_wallCollision(Particle * p){
	if(p->pos.x + p->radius > SCREEN_WIDTH){
		p->pos.x = SCREEN_WIDTH - p->radius;
		p->vel.x*=-1;
	}
	else if(p->pos.x - p->radius < 0){
		p->pos.x = p->radius;
		p->vel.x*=-1;
	}
	
	if(p->pos.y + p->radius > SCREEN_HEIGHT){
		p->pos.y = SCREEN_HEIGHT - p->radius;
		p->vel.y*=-1;
	}
	else if(p->pos.y - p->radius < 0){
		p->pos.y = p->radius;
		p->vel.y*=-1;
	}
}

int particle_pairCollision(Particle * p0, Particle * p1){
	float d = p0->radius + p1->radius;
	if(sqrt((p0->pos.x-p1->pos.x)*(p0->pos.x-p1->pos.x)+(p0->pos.y-p1->pos.y)*(p0->pos.y-p1->pos.y)) < d)return 1;
	return 0;
}
void particle_collision(Particle * p0, Particle * p1){
	//first particle:
	Vector t0;
	vector_subtract(&p0->pos, &p1->pos, &t0);
	Vector t3;
	vector_subtract(&p0->vel, &p1->vel, &t3);
	float k;
	k = vector_dotProduct(&t3,&t0);
	float j;
	j = vector_dotProduct(&t0,&t0);
	float m;
	m = (2*p1->mass)/(p0->mass + p1->mass);
	float u;
	u = (m*k)/j;
	Vector t1;
	vector_multiplyReal(&t0, &t1, u);
	Vector r0;
	vector_subtract(&p0->vel, &t1, &r0);

	//second particle:
	vector_subtract(&p1->pos, &p0->pos, &t0);
	vector_subtract(&p1->vel, &p0->vel, &t3);
	k = vector_dotProduct(&t3,&t0);
	j = vector_dotProduct(&t0,&t0);
	m = (2*p0->mass)/(p0->mass + p1->mass);
	u = (m*k)/j;
	vector_multiplyReal(&t0, &t1, u);
	Vector r1;
	vector_subtract(&p1->vel, &t1, &r1);
	

	vector_multiplyReal(&r0, &p0->vel, 1);
	vector_multiplyReal(&r1, &p1->vel, 1);
	
	//bug handling:
	Vector temp;
	vector_subtract(&p1->pos,&p0->pos,&temp);
	float errd = p1->radius + p0->radius - vector_module(&temp);
	vector_normalize(&temp,&temp);

	Vector temp0, temp1;
	vector_multiplyReal(&temp,&temp0,errd*(p1->mass/(p0->mass+p1->mass)));
	vector_multiplyReal(&temp,&temp1,errd*(p0->mass/(p0->mass+p1->mass)));
	vector_sum(&p1->pos,&temp1,&p1->pos);
	vector_subtract(&p0->pos,&temp0,&p0->pos);

}
float my_rand(float a, float b){
	float r = (float)rand()/RAND_MAX;
	return r*(b-a)+a;
}
int main(int argc, char ** argv) 
{ 
/*
	Vector v0;
	scanf("%f %f", &v0.x, &v0.y);
	Vector v1;
	scanf("%f %f", &v1.x, &v1.y);
	Vector vr;
	vector_projection(&v0,&v1,&vr);
	//vector_rotate(&v0, &vr, PI/2);
	printf("vector op: %f , %f\n",vr.x,vr.y);
	return 0;
*/
	srand(time(0));
	
	//SDL INITIALIZATION
	//
	// retutns zero on success else non-zero 
	if (SDL_Init(SDL_INIT_EVERYTHING) != 0) { 
		printf("error initializing SDL: %s\n", SDL_GetError()); 
	} 
	SDL_Window* win = SDL_CreateWindow("GAME", // creates a window 
									SDL_WINDOWPOS_CENTERED, 
									SDL_WINDOWPOS_CENTERED, 
									SCREEN_WIDTH, SCREEN_HEIGHT, 0); 

	// triggers the program that controls 
	// your graphics hardware and sets flags 
	Uint32 render_flags = SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC; 

	// creates a renderer to render our images 
	SDL_Renderer* rend = SDL_CreateRenderer(win, -1, render_flags); 

	//END OF SDL INITIALIZATION
	
	
	// controls annimation loop 
	int close = 0; 

	// annimation loop 
	int keyboard[1000];
	int keyboardprev[1000];
	int keyboardonce[1000];
	int keyboarddown[1000];
	int keyboardup[1000];
	unsigned int counter;
	int lastTime = SDL_GetTicks();
	memset(keyboard,0,sizeof keyboard);

	Particle p[10000];
	int nparticles=20;
	int i;
	for(i=0;i<nparticles;i++){
		
		p[i].radius = my_rand(10,60);
		
		p[i].color_r = my_rand(0,255);
		p[i].color_g = my_rand(0,255);
		p[i].color_b = my_rand(0,255);
		
		p[i].color_r = 0;
		p[i].color_g = 0;
		p[i].color_b = 0;
		
		//p[i].radius = 20;
		p[i].mass = (4.0/3)*PI*pow(p[i].radius/10.0, 3);
		p[i].pos.x = my_rand(0,SCREEN_WIDTH);
		p[i].pos.y = my_rand(0,SCREEN_HEIGHT);
		float ke = my_rand(0,1000);
		//if(i==0)ke=1000;
		//else ke=0;
		float vel = sqrt((2*ke)/p[i].mass);
		float angle = my_rand(0,2*PI);
		
		p[i].vel.x = vel*cos(angle);
		p[i].vel.y = vel*sin(angle);
	}

	while (!close) {
		SDL_Event event; 
		// Events mangement 
		while (SDL_PollEvent(&event)) {
			switch (event.type) { 

			case SDL_QUIT: 
				// handling of close button 
				close = 1; 
				break; 
			case SDL_KEYDOWN:
				keyboard[event.key.keysym.scancode] = 1;
				keyboarddown[event.key.keysym.scancode] = 1;
				break;
			case SDL_KEYUP:
				keyboard[event.key.keysym.scancode] = 0;
				keyboardup[event.key.keysym.scancode] = 1;
				break;
			}
		} 
		int i;
		for(i=0;i<1000;i++)
		{
			keyboardonce[i] = 0;
			if(!keyboardprev[i] && keyboard[i])keyboardonce[i] = 1;
			keyboardprev[i] = keyboard[i];
		}
		if(keyboardonce[SDL_SCANCODE_RETURN])
		{
		}

		if(keyboardonce[SDL_SCANCODE_K])
		{
		}

		if(keyboardonce[SDL_SCANCODE_L])
		{
		}
		if(keyboard[SDL_SCANCODE_DOWN])
		{
		}
		if(keyboard[SDL_SCANCODE_UP])
		{
		}
		if(keyboard[SDL_SCANCODE_RIGHT])
		{
		}
		if(keyboard[SDL_SCANCODE_LEFT])
		{
		}



		SDL_SetRenderDrawColor( rend, 0, 0, 0, 255 );
		SDL_RenderClear(rend); 
		// clears the screen 

		//Draw circle
		for(i=0;i<nparticles;i++){
			SDL_SetRenderDrawColor( rend, p[i].color_r, p[i].color_g, p[i].color_b, 255 );
			CSDL_DrawCircle(rend, round(p[i].pos.x), round(p[i].pos.y), round(p[i].radius));
			p[i].pos.x+=p[i].vel.x*deltaTime;
			p[i].pos.y+=p[i].vel.y*deltaTime;
			particle_wallCollision(&p[i]);
		}
		int j;
		for(i=0;i<nparticles;i++){
			for(j=i+1;j<nparticles;j++){
				if(particle_pairCollision(&p[i], &p[j])){
					//printf("Colision happening between %d and %d\n",i,j);
					particle_collision(&p[i], &p[j]);
				}
			}
		}
		float kineticenergy=0;
		float totalvel=0;
		float maxvel=vector_module(&p[0].vel);
		float minvel=vector_module(&p[0].vel);
		float maxke=0;
		for(i=0;i<nparticles;i++){
			float vel = vector_module(&p[i].vel);
			float ke = 0.5*p[i].mass*pow(vector_module(&p[i].vel),2);;
			
		//	p[i].color_r = 120*vel;
			//printf("vel = %f\n",vel);
		//	p[i].color_b = 255-p[i].color_r;
			kineticenergy+=0.5*p[i].mass*pow(vector_module(&p[i].vel),2);;
			totalvel += vel;
			if(ke>maxke)maxke=ke;
			if(vel>maxvel)maxvel=vel;
			if(vel<minvel)minvel=vel;
		}
		printf("Total energy = %f, Average speed = %f, maxspeed = %f, minspeed = %f\n",kineticenergy,totalvel/nparticles, maxvel, minvel);
		//Temperature color:
		
		for(i=0;i<nparticles;i++){
			float maxvelBasedOnKE = sqrt(2*kineticenergy/p[i].mass);
			float vel = vector_module(&p[i].vel);
			float ke = 0.5*p[i].mass*pow(vector_module(&p[i].vel),2);;
			p[i].color_r = 255*((ke)/(maxke));
			//printf("vel = %f\n",vel);
			p[i].color_b = 255-p[i].color_r;
			//if(vel == maxvel){
			//	p[i].color_g=200;
			//	p[i].color_r=255;
			//	p[i].color_b=0;
			//}
			//else 
				p[i].color_g=0;
		}
		//particle_wallCollision(&pp);
		//SDL_SetRenderDrawColor( rend, 255, 0, 0, 255 );
		//CSDL_DrawCircle(rend, pp2.pos.x, pp2.pos.y, pp2.radius);
		//pp2.pos.x+=pp2.vel.x;
		//pp2.pos.y+=pp2.vel.y;
		//particle_wallCollision(&pp2);
		//if(particle_pairCollision(&pp, &pp2)){
		//	particle_collision(&pp, &pp2);
			//pp.vel.x=0;
			//pp.vel.y=0;
			//pp2.vel.x=0;
			//pp2.vel.y=0;
		//}
		//textDisplay("Hello", SCREEN_WIDTH/30, SCREEN_WIDTH/2, SCREEN_WIDTH/25, 255, 255, 255, 3);
		
		// triggers the double buffers 
		// for multiple rendering 
		SDL_RenderPresent(rend); 

		// calculates to 60 fps 
		counter++;
		SDL_Delay(1000 / FPS);
	} 

	//SDL CLEANUP
	
	// destroy renderer 
	SDL_DestroyRenderer(rend); 

	// destroy window 
	SDL_DestroyWindow(win);	
	//END OF SDL CLEANUP
	return 0; 
} 

