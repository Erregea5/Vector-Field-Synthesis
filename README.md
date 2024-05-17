# Vector-Field-Synthesis

Vector field synthesis algorithms

Research with Professor Guoning Chen and Graduate Student Nguyen Phan under the Cahsi LREU

Tests use matplotlib to display vector fields.
Vector field data is stored in ply files.
Numpy takes care of the heavy linear algebra calculations


Abstract:
Vector field synthesis from sparse constraints allows for the reconstruction of field-gathered vector flow, efficient and minimally lossy storage of flow data, human designed vector fields for industrial applications, etc. The traditional method of synthesizing vector fields only interpolates the sparse vector constraints and doesn’t consider more physically relevant data such as Jacobian, curl, and divergence, which, depending on the application may be readily available. To address this, we’ve created and tested multiple algorithms with the central focus of finding a balance between performance and accuracy. The methodology of each algorithm revolves around the construction of a linear system via mapping the given constraints to different amounts and types of linear equations. Namely, some methods first augment the sparse field using physical constraints near the known vectors then continue using regular interpolation, others alter the interpolation equations to be more accurate given the physical constraints, and others create an M x N matrix with physical constraints separate from the vector constraints. Indeed, from our results, we observe that including Jacobean constraints improves the accuracy of the algorithm while not having a significant impact on performance. However, we observe that the last type of method suffers from overfitting given too many physical constraints. Additionally, we found that physical constraints too distant from any vector constraints actually end up increasing the error of the system when given too few constraints. These results imply further investigation is in place to optimize the listed hyperparameters.  