#ifndef GPC_COLLISION_INT
#define GPC_COLLISION_INT

namespace Sprog{

	class CollisionIntgl{

		static double Omega11[37*5];
		static double Omega22[37*5];
		static double TStar[37];

	public:
		CollisionIntgl(){}
		~CollisionIntgl(){}
	};
}
#endif


