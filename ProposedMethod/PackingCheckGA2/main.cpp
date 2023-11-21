#include<iostream>
#include<vector>
#include<fstream>
#include<math.h>
#include<sstream>
#include<algorithm>
#include<stdlib.h>
#include<random>
#include<time.h>

using namespace std;

string filename_index = "a1"; //�t�@�C���̖��O
int file_number = 0; //�t�@�C���ԍ�
int total_pieces = 0; //�s�[�X��
double grobal_width = 0;
double grobal_height = 500;
double r1[2][2];
string readfilename = "Jakobs1.CSV";
string readflamename = "flame1.CSV";

//�s�[�X�̃f�[�^�\��
typedef struct pieces {
	int pieceID; //�s�[�XID
	int total_vertex = 0; //���_�̐�
	vector<double> vertex_x; //���_��x���W
	vector<double> vertex_y; //���_��y���W
	double area = 0; //�s�[�X�̖ʐ�
	double graham_area = 0; //�ʖ@�ɂ��ʐ�
	int range_num[2][2] = { 0 }; //[i][j] => [i] : �ŏ��l�E�ő�l [j] : x , y ���W
	bool used = false; //�����`�Ƃ��Ĉ���ꂽ���ǂ���
};
typedef struct flame {
	int pieceID; //�s�[�XID
	int total_vertex = 0; //���_�̐�
	vector<double> vertex_x; //���_��x���W
	vector<double> vertex_y; //���_��y���W
	double area = 0; //�s�[�X�̖ʐ�
	double graham_area = 0; //�ʖ@�ɂ��ʐ�
	int range_num[2][2] = { 0 }; //[i][j] => [i] : �ŏ��l�E�ő�l [j] : x , y ���W
	bool used = false; //�����`�Ƃ��Ĉ���ꂽ���ǂ���
};

//�̂̃f�[�^�\��
typedef struct individual {
	double score = 0;
	vector<double> gen;
};

///////////////////////////////�s�[�X�̊֌W�̊֐�///////////////////////////////////////////////////////////////////////////////////
//�t�@�C���ǂݍ��ݗp
vector<pieces> read_file(vector<pieces> read_piece);
vector<string> split(string& input, char delimiter);
int countLinesInFile(const std::string& filename);

//�s�[�X��x,y���W�̍ő�l�E�ŏ��l���v�Z
vector<pieces> calc_MinMax_vertex(vector<pieces> piecedata);
pieces calc_MinMax_vertex(pieces piecedata);

//�s�[�X�̖ʐς��v�Z(��͑S�Ẵs�[�X�C����1�̃s�[�X)
vector<pieces> calc_piecesArea(vector<pieces> piecedata);
double calc_piecesArea(pieces piecedata);

//�s�[�X�̏�����
vector<pieces> initializePiece(vector<pieces> piecedata);

//�s�[�X�̏o��(��͑S�Ẵs�[�X�C����1�̃s�[�X)
void output_piece_vertex(vector<pieces> output_piece);
void output_piece_vertex(pieces output_piece);

//�s�[�X�̏�������(��͑S�Ẵs�[�X�C����1�̃s�[�X)
void write_pieces_vertex(vector<pieces> write_piece);
void write_pieces_vertex(pieces write_piece);

void removeSpaces(std::string& str);
bool shouldRemove(char c);

pieces Graham_scan(pieces piecedata); //Graham Scan �ɂ��ʖ@
pieces Scan(vector<pieces> p1, pieces p2);
pieces Generate_polygon(vector<int> vertex_id, pieces piecedata); //�ʖ@�̌��ʂ��瑽�p�`�𐶐�
pieces pieceRotation(pieces pieceData, double degree); //�s�[�X����]������
void write(pieces write_piece);
void write_b(vector<pieces> write_piece);

//�ӂ��Ȃ��p���v�Z
double calcDegree(double a1, double a2, double b1, double b2);

//��s�[�X�̕��ƍ�����Ԃ��֐�
vector<double> calcEdge(vector<pieces> data);

bool compareVertexnum(const pieces x, const pieces y); //���_���ɂ��\�[�g�̂��߂̃I�y���[�^
bool compareGraham_Area(const pieces x, const pieces y); //�ʖ@�Ƃ̖ʐύ��ɂ��\�[�g�̂��߂̃I�y���[�^
bool compareScore(const individual x, const individual y);// �X�R�A���Ƀ\�[�g���邽�߂̃I�y���[�^

///////////////////////////��`�I�A���S���Y���֌W/////////////////////////////////////////////////////////////////////////////////////////////////////////////

vector<individual> initializeIndivi(pieces p1, vector<individual> data, int genLength, double minX, double minY, double maxX, double maxY);
vector<individual> Cross(vector<individual> parent, vector<individual> child);
vector<individual> DGA_Cross(vector<individual> parent, vector<individual> child);
vector<individual> Mutation(vector<individual> parent);
bool CrossCheck(pieces p1, pieces p2);
bool InputCheck(pieces p1, pieces p2);
double Evaluate(pieces p1, pieces p2, vector<double> gen);
double Evaluate(vector<pieces> p1, pieces p2, vector<double> gen, vector<double> edge, int sedai);

int main(void) {

	srand((unsigned int)time(NULL));
	bool rectFlag = false; //�s�[�X�Q�������`�Ƃ��Ĉ���ꂽ���ǂ����̃t���O
	for (int aa = 0; aa < 5; aa++) {
		int length = 3;
		int individual_num = 150; //GA�̌̐�
		vector<double> base_edge(2); //��s�[�X�̕��ƍ��������ϐ�
		vector<int> currentPieceID; //�ǂ̃s�[�X���W�c�Ƃ��Ĉ����Ă��邩
		pieces grahamP; // �����̃s�[�X���i�[���邽�߂̉��z�̃s�[�X�ϐ��i�ʕ�ɗ��p���邽�߁j
		vector<individual> parents(individual_num);
		vector<individual> children(individual_num);
		int total_island = 8;
		int generation = 200; //GA�̐��㐔
		pieces pair2;
		int total_pieces = countLinesInFile(readfilename);
		vector<pieces> fixedPiece(total_pieces); //�s�[�X�f�[�^�Z�b�g�i�ŏI�I�ȉ��ƂȂ�悤�Ɉʒu�𓮂����Ă����j
		vector<pieces> flame(1); //�g�̓ǂݍ���
		vector<pieces> pairPiece, pair1;
		fixedPiece = initializePiece(fixedPiece);
		//write_pieces_vertex(fixedPiece);
		double areaSum = 0;
		for (int i = 0; i < (int)fixedPiece.size(); i++) areaSum += fixedPiece.at(i).area;
		//�����l�̌���
		pairPiece = fixedPiece;
		cout << areaSum << endl;
		sort(pairPiece.begin(), pairPiece.end(), compareGraham_Area);
		pair1.push_back(pairPiece.at(0)), currentPieceID.push_back(pairPiece.at(0).pieceID);
		r1[0][0] = pair1.at(0).vertex_x.at(pair1.at(0).range_num[0][0]);
		r1[0][1] = pair1.at(0).vertex_y.at(pair1.at(0).range_num[0][1]);
		r1[1][0] = pair1.at(0).vertex_x.at(pair1.at(0).range_num[1][0]);
		r1[1][1] = pair1.at(0).vertex_y.at(pair1.at(0).range_num[1][1]);
		grobal_width = abs(r1[0][0] - r1[1][0]);
		//�����`�̐���
		for (int s = 0; (int)s < total_pieces - 1; s++) {
			pair2 = pairPiece.at(total_pieces - s - 1);
			if (rectFlag) vector<int> currentPieceID;
			parents = initializeIndivi(pair2, parents, length, r1[0][0], r1[0][1], r1[1][0], r1[1][1]);
			children = initializeIndivi(pair2, parents, length, r1[0][0], r1[0][1], r1[1][0], r1[1][1]);
			for (int j = 0; j < (int)parents.size(); j++) parents.at(j).score = Evaluate(pair1, pair2, parents.at(j).gen, base_edge, 0);
			for (int j = 0; j < total_island; j++) {
				if (j == 0) sort(parents.begin(), parents.begin() + total_island, compareScore);
				else if (j == total_island - 1) sort(parents.begin() + j * total_island, parents.end(), compareScore);
				else sort(parents.begin() + j * total_island, parents.begin() + (j + 1) * total_island, compareScore);
			}

			base_edge = calcEdge(pair1);
			//��`�I�A���S���Y���K�p
			for (int i = 0; i < generation; i++) {
				parents = DGA_Cross(parents, children);
				parents = Mutation(parents);
				for (int j = 0; j < (int)parents.size(); j++) parents.at(j).score = Evaluate(pair1, pair2, parents.at(j).gen, base_edge, i);
				sort(parents.begin(), parents.end(), compareScore);
			}
			//��`�q�������ƂɑI�������s�[�X�̈ړ��Ɖ�]���s��
			for (int i = 0; i < fixedPiece.at(pair2.pieceID).total_vertex; i++) {
				fixedPiece.at(pair2.pieceID).vertex_x.at(i) += parents.at(0).gen.at(0);
				fixedPiece.at(pair2.pieceID).vertex_y.at(i) += parents.at(0).gen.at(1);
			}
			fixedPiece.at(pair2.pieceID) = pieceRotation(fixedPiece.at(pair2.pieceID), parents.at(0).gen.at(2));
			fixedPiece.at(pair2.pieceID) = calc_MinMax_vertex(fixedPiece.at(pair2.pieceID));
			r1[0][0] = min(r1[0][0], fixedPiece.at(pair2.pieceID).vertex_x.at(fixedPiece.at(pair2.pieceID).range_num[0][0]));
			r1[0][1] = min(r1[0][1], fixedPiece.at(pair2.pieceID).vertex_y.at(fixedPiece.at(pair2.pieceID).range_num[0][1]));
			r1[1][0] = max(r1[1][0], fixedPiece.at(pair2.pieceID).vertex_x.at(fixedPiece.at(pair2.pieceID).range_num[1][0]));
			r1[1][1] = max(r1[1][1], fixedPiece.at(pair2.pieceID).vertex_y.at(fixedPiece.at(pair2.pieceID).range_num[1][1]));
			grobal_width = abs(r1[0][0] - r1[1][0]);
			pair1.push_back(fixedPiece.at(pair2.pieceID));
			write_b(pair1);
			grahamP = Scan(pair1, pair2);
			currentPieceID.push_back(grahamP.pieceID);

			write(grahamP);
			int point = 1;
			for (int i = 0; i < (int)flame.size(); i++) if (InputCheck(flame.at(i), grahamP)) point = 0;

			cout << point << endl;
			if (point == 0) {
				cout << "OK" << endl;
			}
		}
		write_pieces_vertex(fixedPiece);
		cout << "aaaaaa" << endl;
	}
}

/////////////////////////////////////////////////////�t�@�C���̓ǂݍ��݊֌W/////////////////////////////////////////////////////////////////////////////////////////
//�s���𐔂��邽�߂����̂��
int countLinesInFile(const std::string& filename) {
	std::ifstream file(filename);
	if (!file.is_open()) {
		std::cerr << "Error opening file: " << filename << std::endl;
		return -1;
	}

	int lineCount = 0;
	std::string line;
	while (std::getline(file, line)) {
		lineCount++;
	}

	return lineCount;
}
//�t�@�C���̐��l����
vector<pieces> read_file(vector<pieces> read_piece) {
	ifstream ifs(readfilename);
	string line;
	int count = 0, i = 0;
	while (getline(ifs, line)) {
		vector<string> strvec = split(line, ',');
		//if (strvec.size() % 2 != 0 && strvec.size() >= 3) {
			//strvec.erase(strvec.end()-1);
		//}
		for (auto it = strvec.begin(); it != strvec.end();) {
			if ((*it).empty()) {
				it = strvec.erase(it);
			}
			else {
				++it;
			}
		}
		//for (std::string& str : strvec) {
			//removeSpaces(str);
		//}

		for (i = 0; i < (int)strvec.size(); i++) {
			if (i % 2 == 0) read_piece.at(count).vertex_x.push_back(stoi(strvec.at(i)));
			else read_piece.at(count).vertex_y.push_back(stoi(strvec.at(i)));
		}
		read_piece.at(count).total_vertex = (i + 1) / 2;
		read_piece.at(count).pieceID = count;
		count++;

	}
	total_pieces = count;
	cout << total_pieces << endl;
	return read_piece;
}
//�t�@�C����ǂݍ���
vector<string> split(string& input, char delimiter) {
	istringstream stream(input);
	string field;
	vector<string> result;
	while (getline(stream, field, delimiter)) {
		result.push_back(field);
	}
	return result;
}

//void removeSpaces(std::string& str) {
	//str.erase(std::remove_if(str.begin(), str.end(), shouldRemove), str.end());
//}
//bool shouldRemove(char c) {
	//return !std::isspace(c) && !std::isdigit(c);
//}
///////////////////////////////////�s�[�X��x,y���W�̍ő�E�ŏ����_�ԍ����v�Z///////////////////////////////////////////////////////////////////////////////////////
vector<pieces> calc_MinMax_vertex(vector<pieces> piecedata) {
	vector<vector<double>> range(2, vector<double>(2)); //range[i][j] : i �͍ő�E�ŏ��Cj ��x,y���W (i = 0�̂Ƃ��ŏ�, 1�̂Ƃ��ő�)
	for (int i = 0; i < (int)piecedata.size(); i++) {
		range.at(0).at(0) = piecedata.at(i).vertex_x.at(0);
		range.at(1).at(0) = piecedata.at(i).vertex_x.at(0);
		range.at(0).at(1) = piecedata.at(i).vertex_y.at(0);
		range.at(1).at(1) = piecedata.at(i).vertex_y.at(0);
		for (int j = 0; j < piecedata.at(i).total_vertex; j++) {
			if (range.at(0).at(0) >= piecedata.at(i).vertex_x.at(j)) {
				range.at(0).at(0) = piecedata.at(i).vertex_x.at(j);
				piecedata.at(i).range_num[0][0] = j;
			}
			if (range.at(1).at(0) <= piecedata.at(i).vertex_x.at(j)) {
				range.at(1).at(0) = piecedata.at(i).vertex_x.at(j);
				piecedata.at(i).range_num[1][0] = j;
			}
			if (range.at(0).at(1) >= piecedata.at(i).vertex_y.at(j)) {
				range.at(0).at(1) = piecedata.at(i).vertex_y.at(j);
				piecedata.at(i).range_num[0][1] = j;
			}
			if (range.at(1).at(1) <= piecedata.at(i).vertex_y.at(j)) {
				range.at(1).at(1) = piecedata.at(i).vertex_y.at(j);
				piecedata.at(i).range_num[1][1] = j;
			}
		}
	}
	return piecedata;
}

pieces calc_MinMax_vertex(pieces piecedata) {
	vector<vector<double>> range(2, vector<double>(2)); //range[i][j] : i �͍ő�E�ŏ��Cj ��x,y���W (i = 0�̂Ƃ��ŏ�, 1�̂Ƃ��ő�)
	range.at(0).at(0) = piecedata.vertex_x.at(0);
	range.at(1).at(0) = piecedata.vertex_x.at(0);
	range.at(0).at(1) = piecedata.vertex_y.at(0);
	range.at(1).at(1) = piecedata.vertex_y.at(0);
	for (int j = 0; j < piecedata.total_vertex; j++) {
		if (range.at(0).at(0) >= piecedata.vertex_x.at(j)) {
			range.at(0).at(0) = piecedata.vertex_x.at(j);
			piecedata.range_num[0][0] = j;
		}
		if (range.at(1).at(0) <= piecedata.vertex_x.at(j)) {
			range.at(1).at(0) = piecedata.vertex_x.at(j);
			piecedata.range_num[1][0] = j;
		}
		if (range.at(0).at(1) >= piecedata.vertex_y.at(j)) {
			range.at(0).at(1) = piecedata.vertex_y.at(j);
			piecedata.range_num[0][1] = j;
		}
		if (range.at(1).at(1) <= piecedata.vertex_y.at(j)) {
			range.at(1).at(1) = piecedata.vertex_y.at(j);
			piecedata.range_num[1][1] = j;
		}
	}
	return piecedata;
}


//�s�[�X�̖ʐς��v�Z////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<pieces> calc_piecesArea(vector<pieces > piecedata) {
	for (int i = 0; i < (int)piecedata.size(); i++) {
		for (int j = 0; j < piecedata.at(i).total_vertex; j++) {
			if (j == piecedata.at(i).total_vertex - 1) {
				piecedata.at(i).area += piecedata.at(i).vertex_x.at(j) * piecedata.at(i).vertex_y.at(0) - piecedata.at(i).vertex_y.at(j) * piecedata.at(i).vertex_x.at(0);
			}
			else {
				piecedata.at(i).area += piecedata.at(i).vertex_x.at(j) * piecedata.at(i).vertex_y.at(j + 1) - piecedata.at(i).vertex_y.at(j) * piecedata.at(i).vertex_x.at(j + 1);
			}
		}
		piecedata.at(i).area = 0.5 * abs(piecedata.at(i).area);
	}
	return piecedata;
}

double calc_piecesArea(pieces piecedata) {
	for (int j = 0; j < piecedata.total_vertex; j++) {
		if (j == piecedata.total_vertex - 1) {
			piecedata.area += piecedata.vertex_x.at(j) * piecedata.vertex_y.at(0) - piecedata.vertex_y.at(j) * piecedata.vertex_x.at(0);
		}
		else {
			piecedata.area += piecedata.vertex_x.at(j) * piecedata.vertex_y.at(j + 1) - piecedata.vertex_y.at(j) * piecedata.vertex_x.at(j + 1);
		}
	}
	piecedata.area = 0.5 * abs(piecedata.area);
	return piecedata.area;
}

////////////////////////////////////////////////�s�[�X�̏��������s��/////////////////////////////////////////////////////////////////////////////////////////////////
vector<pieces> initializePiece(vector<pieces> piecedata) {
	piecedata = read_file(piecedata);
	piecedata = calc_MinMax_vertex(piecedata);
	piecedata = calc_MinMax_vertex(piecedata);
	piecedata = calc_piecesArea(piecedata);
	for (int i = 0; i < (int)piecedata.size(); i++) {
		piecedata.at(i).graham_area = calc_piecesArea(Graham_scan(piecedata.at(i))) - piecedata.at(i).area;
	}
	return piecedata;
}

////////////////////////////////////////////////�e�s�[�X�̒��_���W��\��//////////////////////////////////////////////////////////////////////////////////////////////
void output_piece_vertex(vector<pieces> output_piece) {
	for (int i = 0; i < (int)output_piece.size(); i++) {
		cout << endl << "Piece" << i << endl;
		for (int j = 0; j < output_piece.at(i).total_vertex; j++) cout << "(" << output_piece.at(i).vertex_x.at(j) << " , " << output_piece.at(i).vertex_y.at(j) << ")" << endl;
	}
}

void output_piece_vertex(pieces output_piece) {
	cout << endl << "Piece" << endl;
	for (int j = 0; j < output_piece.total_vertex; j++) cout << "(" << output_piece.vertex_x.at(j) << " , " << output_piece.vertex_y.at(j) << ")" << endl;
}


////////////////////////////////�e�s�[�X�̒��_���W����������//////////////////////////////////////////////////////////////////////////////////////////////////////////
void write_pieces_vertex(vector<pieces> write_piece) {
	string filename = filename_index + to_string(file_number) + ".CSV";
	ofstream ofs(filename);
	for (int i = 0; i < (int)write_piece.size(); i++) {
		for (int j = 0; (int)j < write_piece.at(i).total_vertex; j++) {
			ofs << write_piece.at(i).vertex_x.at(j) << "," << write_piece.at(i).vertex_y.at(j) << "," << endl;
		}
		ofs << write_piece.at(i).vertex_x.at(0) << "," << write_piece.at(i).vertex_y.at(0) << "," << endl << endl;
	}
	file_number++;
}

void write(pieces write_piece) {
	string filename = "after" + to_string(file_number) + ".CSV";
	ofstream ofs(filename);
	for (int j = 0; (int)j < write_piece.total_vertex; j++) {
		ofs << write_piece.vertex_x.at(j) << "," << write_piece.vertex_y.at(j) << "," << endl;
	}
	ofs << write_piece.vertex_x.at(0) << "," << write_piece.vertex_y.at(0) << "," << endl << endl;
	file_number++;
}
void write_b(vector<pieces> write_piece) {
	string filename = filename_index + to_string(file_number) + ".CSV";
	ofstream ofs(filename);
	for (int i = 0; i < (int)write_piece.size(); i++) {
		for (int j = 0; (int)j < write_piece.at(i).total_vertex; j++) {
			ofs << write_piece.at(i).vertex_x.at(j) << "," << write_piece.at(i).vertex_y.at(j) << "," << endl;
		}
		ofs << write_piece.at(i).vertex_x.at(0) << "," << write_piece.at(i).vertex_y.at(0) << "," << endl << endl;
	}
	file_number++;
}


void write_pieces_vertex(pieces write_piece) {
	string filename = filename_index + to_string(file_number) + ".CSV";
	ofstream ofs(filename);
	for (int j = 0; (int)j < write_piece.total_vertex; j++) {
		ofs << write_piece.vertex_x.at(j) << "," << write_piece.vertex_y.at(j) << "," << endl;
	}
	ofs << write_piece.vertex_x.at(0) << "," << write_piece.vertex_y.at(0) << "," << endl << endl;
	file_number++;
}

///////////////////////////////Graham_scan(�ʖ@)//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
pieces Graham_scan(pieces piecedata) {
	vector<vector<double>> point(piecedata.total_vertex - 1, vector<double>(3)); //���p�`�̒��_�f�[�^(�s�[�X�ԍ��̒���{���ρCx���W�CID})
	vector<double> a(2), b(2); //�x�N�g��a,b�𐶐�
	double abs_a, abs_b; //�x�N�g���̐�Βl
	//y���W���ŏ��ƂȂ钸�_�̍��W�𐶐�
	double x0 = piecedata.vertex_x.at(piecedata.range_num[0][1]);
	double y0 = piecedata.vertex_y.at(piecedata.range_num[0][1]);
	//���ς��v�Z
	for (int i = 0; i < (int)piecedata.total_vertex; i++) {
		int j = i;
		if (i == piecedata.range_num[0][1]) continue;
		else if (i > piecedata.range_num[0][1]) j--;
		point.at(j).at(1) = piecedata.vertex_x.at(i);
		point.at(j).at(2) = (double)i;
		a.at(0) = piecedata.vertex_x.at(i) - x0;
		a.at(1) = piecedata.vertex_y.at(i) - y0;
		b.at(0) = 1;
		b.at(1) = 0;
		abs_a = sqrt(pow(a.at(0), 2) + pow(a.at(1), 2));
		abs_b = sqrt(pow(b.at(0), 2) + pow(b.at(1), 2));
		if (abs_a * abs_b != 0) point.at(j).at(0) = (a.at(0) * b.at(0) + a.at(1) * b.at(1)) / (abs_a * abs_b);
		else cout << "Error : Abs a or b is 0" << endl;
	}
	sort(point.begin(), point.end());
	reverse(point.begin(), point.end());
	vector<int> stack_vertex;
	//�X�^�b�N�̏�����
	stack_vertex.push_back(piecedata.range_num[0][1]);
	stack_vertex.push_back((int)point.at(0).at(2));
	stack_vertex.push_back((int)point.at(1).at(2));
	int max_num = 2;
	while (1) {
		if (stack_vertex.at((int)stack_vertex.size() - 1) == point.at((int)point.size() - 1).at(2)) break;
		double Px = piecedata.vertex_x.at(stack_vertex.at((int)stack_vertex.size() - 2)), Py = piecedata.vertex_y.at(stack_vertex.at((int)stack_vertex.size() - 2));
		double Cx = piecedata.vertex_x.at(stack_vertex.at((int)stack_vertex.size() - 1)), Cy = piecedata.vertex_y.at(stack_vertex.at((int)stack_vertex.size() - 1));
		double Qx = piecedata.vertex_x.at((int)point.at(max_num).at(2)), Qy = piecedata.vertex_y.at((int)point.at(max_num).at(2));
		if (((Px - Cx) * (Qy - Cy) - (Py - Cy) * (Qx - Cx)) < 0) {
			stack_vertex.push_back((int)point.at(max_num).at(2));
			max_num++;
		}
		else {
			stack_vertex.pop_back();
			if ((int)stack_vertex.size() <= 2) {
				stack_vertex.push_back((int)point.at(max_num).at(2));
				max_num++;
			}
		}
	}
	return Generate_polygon(stack_vertex, piecedata);
}

pieces Scan(vector<pieces> p1, pieces p2) {
	pieces grahamP; // �����̃s�[�X���i�[���邽�߂̉��z�̃s�[�X�ϐ��i�ʕ�ɗ��p���邽�߁j
	double total_area = p2.area, width, height;
	int add_total_vertex = 0; //���₵�����_�̃J�E���^

	grahamP = p2;
	for (int i = 0; i < (int)p1.size(); i++) {
		for (int j = 0; j < p1.at(i).total_vertex; j++) {
			grahamP.vertex_x.push_back(p1.at(i).vertex_x.at(j));
			grahamP.vertex_y.push_back(p1.at(i).vertex_y.at(j));
			add_total_vertex++;
		}
	}
	grahamP.total_vertex += add_total_vertex;
	grahamP = calc_MinMax_vertex(grahamP);
	return Graham_scan(grahamP);
}
//�ʖ@�̌��ʂ��瑽�p�`�𐶐�
pieces Generate_polygon(vector<int> vertex_id, pieces  piecedata) {
	pieces polygon;
	polygon.total_vertex = (int)vertex_id.size();
	for (int j = 0; j < (int)vertex_id.size(); j++) {
		polygon.vertex_x.push_back(piecedata.vertex_x.at(vertex_id.at(j)));
		polygon.vertex_y.push_back(piecedata.vertex_y.at(vertex_id.at(j)));
	}
	return polygon;
}

//�s�[�X���p�xdegree�ŉ�]////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
pieces pieceRotation(pieces pieceData, double degree) {
	const double PI = 3.1415926535893238;
	double temp_x, temp_y;
	degree = degree * PI / 180;
	for (int i = 0; i < pieceData.total_vertex; i++) {
		temp_x = pieceData.vertex_x.at(i), temp_y = pieceData.vertex_y.at(i);
		pieceData.vertex_x.at(i) = temp_x * cos(degree) - temp_y * sin(degree);
		pieceData.vertex_y.at(i) = temp_x * sin(degree) + temp_y * cos(degree);
	}
	pieceData = calc_MinMax_vertex(pieceData);
	return pieceData;
}

//�ӓ��m���Ȃ��p���v�Z////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double calcDegree(double a1, double a2, double b1, double b2) {
	double abs_a = sqrt(pow(a1, 2) + pow(a2, 2));
	double abs_b = sqrt(pow(b1, 2) + pow(b2, 2));
	double cosDeg = (a1 * b1 + a2 * b2) / (abs_a * abs_b);
	return acos(cosDeg) * 180 / 3.141592653589793238;
}


//��s�[�X�̕��ƍ������v�Z////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<double> calcEdge(vector<pieces> data) {
	vector<double> width(2); //0�����͍ŏ���x���W�C1�����͍ő��x���W
	vector<double> height(2); //0�����͍ŏ�y���W�C1�����͍ő��y���W
	vector<double> edge(2); //��s�[�X�̕��ƍ�������
	//�ϐ��̏�����
	width.at(0) = data.at(0).vertex_x.at(data.at(0).range_num[0][0]);
	width.at(1) = data.at(0).vertex_x.at(data.at(0).range_num[1][0]);
	height.at(0) = data.at(0).vertex_y.at(data.at(0).range_num[0][1]);
	height.at(1) = data.at(0).vertex_y.at(data.at(0).range_num[1][1]);
	for (int i = 0; i < (int)data.size(); i++) {
		width.at(0) = min(width.at(0), data.at(i).vertex_x.at(data.at(i).range_num[0][0]));
		width.at(1) = max(width.at(1), data.at(i).vertex_x.at(data.at(i).range_num[1][0]));
		height.at(0) = min(height.at(0), data.at(i).vertex_y.at(data.at(i).range_num[0][1]));
		height.at(1) = min(height.at(1), data.at(i).vertex_y.at(data.at(i).range_num[1][1]));
	}
	edge.at(0) = width.at(1) - width.at(0);
	edge.at(1) = height.at(1) - height.at(0);
	return edge;
}

//�\�[�g�p�̃I�y���[�^
bool compareVertexnum(const pieces x, const pieces y) {
	return x.total_vertex > y.total_vertex;
}

bool compareGraham_Area(const pieces x, const pieces y) {
	return x.graham_area > y.graham_area;
}

bool compareScore(const individual x, const individual y) {
	return x.score > y.score;
}


////////////////////////�̂ƈ�`�q�̏�����///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<individual> initializeIndivi(pieces p1, vector<individual> data, int genLength, double minX, double minY, double maxX, double maxY) {
	double piece_edge[2][2];
	int total_island = 8;
	int count = 0, sub_population = (int)data.size() / total_island;
	piece_edge[0][0] = p1.vertex_x.at(p1.range_num[0][0]); //�I���s�[�X�̍ŏ�x���W
	piece_edge[0][1] = p1.vertex_y.at(p1.range_num[0][1]); //�I���s�[�X�̍ŏ�y���W
	piece_edge[1][0] = p1.vertex_x.at(p1.range_num[1][0]); //�I���s�[�X�̍ő�x���W
	piece_edge[1][1] = p1.vertex_y.at(p1.range_num[1][1]); //�I���s�[�X�̍ő�y���W

	random_device seed_gen;
	mt19937 engine(seed_gen());
	uniform_real_distribution<> norm(-1, 1);
	/*uniform_real_distribution<> initialize_X1(minX - 2 * piece_edge[1][0], minX - piece_edge[1][0]);
	uniform_real_distribution<> initialize_X2(minX, maxX);
	uniform_real_distribution<> initialize_X3(maxX + piece_edge[0][0], maxX + 2 * piece_edge[0][0]);
	uniform_real_distribution<> initialize_Y1(minY - 2 * piece_edge[1][1], minY - piece_edge[1][1]);
	uniform_real_distribution<> initialize_Y2(minY, maxY);
	uniform_real_distribution<> initialize_Y3(maxY + piece_edge[0][1], maxY + 2 * piece_edge[0][1]);*/
	double randnum;
	for (int i = 0; i < (int)data.size(); i++) {
		data.at(i).gen.resize(genLength);
		if (i % sub_population == 0) count++;
		for (int j = 0; j < genLength; j++) {
			randnum = norm(engine);
			if (j == 0) {
				if (count == 1 || count == 4 || count == 6) {
					if (randnum < 0) data.at(i).gen.at(j) = abs(randnum) * (minX - 1.5 * piece_edge[1][0]);
					else data.at(i).gen.at(j) = randnum * (minX - piece_edge[1][0]);
				}
				else if (count == 2 || count == 7) {
					if (randnum < 0) data.at(i).gen.at(j) = abs(randnum) * minX;
					else data.at(i).gen.at(j) = randnum * maxX;
				}
				else {
					if (randnum < 0) data.at(i).gen.at(j) = abs(randnum) * (maxX + piece_edge[0][0]);
					else data.at(i).gen.at(j) = randnum * (maxX + 1.5 * piece_edge[0][0]);
				}
			}
			else if (j == 1) {
				if (count == 1 || count == 2 || count == 3) {
					if (randnum < 0) data.at(i).gen.at(j) = abs(randnum) * (minY - 1.5 * piece_edge[1][1]);
					else data.at(i).gen.at(j) = randnum * (minY - piece_edge[1][1]);
				}
				else if (count == 4 || count == 5) {
					if (randnum < 0) data.at(i).gen.at(j) = abs(randnum) * minY;
					else data.at(i).gen.at(j) = randnum * maxY;
				}
				else {
					if (randnum < 0) data.at(i).gen.at(j) = abs(randnum) * (maxY + piece_edge[0][1]);
					else data.at(i).gen.at(j) = abs(randnum) * (maxY + 1.5 * piece_edge[0][1]);
				}
			}
			else  data.at(i).gen.at(j) = ((double)(rand() % 36) * (double)(j / 10 + 1)); //(double)(rand() % 365);
		}
	}
	return data;
}

//////////////////////BLX-���ɂ�����///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<individual> Cross(vector<individual> parent, vector<individual> child) {
	random_device seed_gen;
	mt19937 engine1(seed_gen());
	uniform_real_distribution<> Generate(-1, 1);
	double alpha = 0.3, value_x, value_y, normDeg;
	vector<vector <double>>  range(2, vector<double>((int)parent.at(0).gen.size())); //1�����ڂ�Min,Max�C2�����ڂ�x,y������
	vector<double> distance((int)parent.at(0).gen.size());
	vector<double>  total_score((int)parent.size(), 0);
	for (int i = 0; i < (int)child.size() / 10; i++) child.at(i) = parent.at(i);
	total_score.at(0) = parent.at(0).score;
	for (int i = 1; i < (int)parent.size(); i++) total_score.at(i) = total_score.at(i - 1) + parent.at(i).score;
	if ((int)total_score.at((int)total_score.size() - 1) > 0) {
		//����
		for (int i = (int)child.size() / 10; i < (int)child.size(); i++) {
			//���[���b�g�I��
			int num1 = 0, num2 = 0;
			double percent1 = rand() % (int)total_score.at((int)total_score.size() - 1);
			double percent2 = rand() % (int)total_score.at((int)total_score.size() - 1);
			while (percent1 >= total_score.at(num1)) num1++;
			while (percent2 >= total_score.at(num2)) num2++;
			//�u�����h����
			for (int j = 0; j < 2; j++) {
				for (int k = 0; k < (int)range.at(j).size(); k++) {
					if (j == 0) {
						distance.at(k) = abs(parent.at(num1).gen.at(k) - parent.at(num2).gen.at(k));
						range.at(j).at(k) = min(parent.at(num1).gen.at(k), parent.at(num2).gen.at(k)) - distance.at(k) * alpha;
					}
					else  range.at(j).at(k) = max(parent.at(num1).gen.at(k), parent.at(num2).gen.at(k)) - distance.at(k) * alpha;
				}
			}
			value_x = Generate(engine1);
			normDeg = abs(parent.at(num1).gen.at(2) - parent.at(num2).gen.at(2));
			child.at(i).gen.at(2) = (min(parent.at(num1).gen.at(2), parent.at(num2).gen.at(2)) + normDeg * 0.5) + Generate(engine1) * normDeg;
			if (child.at(i).gen.at(2) >= 360) child.at(i).gen.at(2) -= 360;
			else if (child.at(i).gen.at(2) < 0) child.at(i).gen.at(2) += 360;
			if (value_x < 0) value_x = value_x * range.at(0).at(0);
			else value_x = value_x * range.at(1).at(0);
			value_y = Generate(engine1);
			if (value_y < 0) value_y = value_y * range.at(0).at(1);
			else value_y = value_y * range.at(1).at(1);
			for (int j = 0; j < (int)child.at(i).gen.size() - 1; j++) {
				if (j == 0) child.at(i).gen.at(j) = value_x;
				else if (j == 1) child.at(i).gen.at(j) = value_y;
			}
		}
	}
	for (int i = 0; i < (int)parent.size(); i++) parent.at(i) = child.at(i);
	return parent;
}

vector<individual> DGA_Cross(vector<individual> parent, vector<individual> child) {
	int island_num[2] = { 0, 8 }; //���̔ԍ�
	int total_island = 8; //���̑���
	random_device seed_gen;
	mt19937 engine1(seed_gen());
	uniform_real_distribution<> Generate(-1, 1);
	double alpha = 0.3, value_x, value_y, normDeg;
	vector<vector <double>>  range(2, vector<double>((int)parent.at(0).gen.size())); //1�����ڂ�Min,Max�C2�����ڂ�x,y������
	vector<double> distance((int)parent.at(0).gen.size()); //�u�����h�����̕�
	vector<double> total_score((int)parent.size() / total_island, 0); //���[���b�g�I���̂��߂ɗݐϘa���i�[����z��
	//DGA�̃T�u��W�c�̌���
	for (int isnum = 0; isnum < 8; isnum++) {
		for (int i = 0; i < (int)child.size() / 10 / total_island; i++) child.at(i) = parent.at(i + isnum * 8);
		total_score.at(0) = parent.at(isnum * total_island).score;
		for (int i = 1; i < (int)parent.size() / total_island; i++) total_score.at(i) = total_score.at(i - 1) + parent.at(i).score;
		if ((int)total_score.at((int)total_score.size() - 1) > 0) {
			for (int i = (int)child.size() / 10 / total_island; i < (int)child.size() / total_island; i++) {
				int num1 = 0, num2 = 0;
				double percent1 = rand() % (int)total_score.at((int)total_score.size() - 1);
				double percent2 = rand() % (int)total_score.at((int)total_score.size() - 1);
				while (percent1 >= total_score.at(num1)) num1++;
				while (percent2 >= total_score.at(num2)) num2++;
				//�u�����h����
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < (int)range.at(j).size(); k++) {
						if (j == 0) {
							distance.at(k) = abs(parent.at(num1 + isnum * total_island).gen.at(k) - parent.at(num2 + isnum * total_island).gen.at(k));
							range.at(j).at(k) = min(parent.at(num1 + isnum * total_island).gen.at(k), parent.at(num2 + isnum * total_island).gen.at(k)) - distance.at(k) * alpha;
						}
						else  range.at(j).at(k) = max(parent.at(num1 + isnum * total_island).gen.at(k), parent.at(num2 + isnum * total_island).gen.at(k)) - distance.at(k) * alpha;
					}
				}
				value_x = Generate(engine1);
				normDeg = abs(parent.at(num1 + isnum * total_island).gen.at(2) - parent.at(num2 + isnum * total_island).gen.at(2));
				child.at(i + isnum * total_island).gen.at(2) = (min(parent.at(num1 + isnum * total_island).gen.at(2), parent.at(num2 + isnum * total_island).gen.at(2)) + normDeg * 0.5) + Generate(engine1) * normDeg;
				if (child.at(i + isnum * total_island).gen.at(2) >= 360) child.at(i + isnum * total_island).gen.at(2) -= 360;
				else if (child.at(i + isnum * total_island).gen.at(2) < 0) child.at(i + isnum * total_island).gen.at(2) += 360;
				if (value_x < 0) value_x = value_x * range.at(0).at(0);
				else value_x = value_x * range.at(1).at(0);
				value_y = Generate(engine1);
				if (value_y < 0) value_y = value_y * range.at(0).at(1);
				else value_y = value_y * range.at(1).at(1);
				for (int j = 0; j < (int)child.at(i + isnum * total_island).gen.size() - 1; j++) {
					if (j == 0) child.at(i + isnum * total_island).gen.at(j) = value_x;
					else if (j == 1) child.at(i + isnum * total_island).gen.at(j) = value_y;
				}
			}
		}
		for (int i = 0; i < (int)parent.size() / total_island; i++) parent.at(i + isnum * total_island) = child.at(i + isnum * total_island);
	}
	return parent;
}

///////////�]���֐�////////////////////////////////////////////////////////////////////////////////////////////////////
double Evaluate(pieces p1, pieces p2, vector<double> gen) {
	double score, total_area = p1.area + p2.area;
	int add_total_vertex = 0;
	for (int i = 0; i < (int)p2.total_vertex; i++) {
		p2.vertex_x.at(i) += gen.at(0);
		p2.vertex_y.at(i) += gen.at(1);
	}
	p2 = pieceRotation(p2, gen.at(2));
	if (CrossCheck(p1, p2)) return 0;//���㐔���㔼�ɂȂ�قǌ������O���͂�邭�iDGA�j�ŏd�Ȃ�ɂ���
	for (int i = 0; i < (int)p2.total_vertex; i++) {
		p1.vertex_x.push_back(p2.vertex_x.at(i));
		p1.vertex_y.push_back(p2.vertex_y.at(i));
		add_total_vertex++;
	}
	p1.total_vertex += add_total_vertex;
	p1 = calc_MinMax_vertex(p1);
	score = calc_piecesArea(Graham_scan(p1));
	score = pow(total_area / score, 3);
	return score * 100;//3��100�͑��삵�Ď����̂����l�Ɍ��߂�
}

double Evaluate(vector<pieces> p1, pieces p2, vector<double> gen, vector<double> edge, int sedai) {
	double score, score2, score3, total_area = p2.area, width, height;
	int add_total_vertex = 0; //���₵�����_�̃J�E���^
	pieces grahamP; // �����̃s�[�X���i�[���邽�߂̉��z�̃s�[�X�ϐ��i�ʕ�ɗ��p���邽�߁j
	for (int i = 0; i < (int)p1.size(); i++) total_area += p1.at(i).area;
	for (int i = 0; i < (int)p2.total_vertex; i++) {
		p2.vertex_x.at(i) += gen.at(0);
		p2.vertex_y.at(i) += gen.at(1);
	}
	p2 = pieceRotation(p2, gen.at(2));
	p2 = calc_MinMax_vertex(p2);
	//width = abs(p2.vertex_x.at(p2.range_num[0][0]) - p2.vertex_x.at(p2.range_num[1][0]));
	grahamP = p2;
	for (int i = 0; i < (int)p1.size(); i++) if (CrossCheck(p1.at(i), p2)) return 0;
	for (int i = 0; i < (int)p1.size(); i++) if (InputCheck(p1.at(i), p2)) return 0;
	for (int i = 0; i < (int)p1.size(); i++) {
		for (int j = 0; j < p1.at(i).total_vertex; j++) {
			grahamP.vertex_x.push_back(p1.at(i).vertex_x.at(j));
			grahamP.vertex_y.push_back(p1.at(i).vertex_y.at(j));
			add_total_vertex++;
		}
	}
	grahamP.total_vertex += add_total_vertex;
	grahamP = calc_MinMax_vertex(grahamP);
	width = abs(grahamP.vertex_x.at(grahamP.range_num[0][0]) - grahamP.vertex_x.at(grahamP.range_num[1][0])); //�g��������
	height = abs(grahamP.vertex_y.at(grahamP.range_num[0][1]) - grahamP.vertex_y.at(grahamP.range_num[1][1])); //�g��������
	if (height > grobal_height) return 0;
	score = total_area / calc_piecesArea(Graham_scan(grahamP));
	score2 = width * height;
	score2 = total_area / score2;
	score3 = grobal_width / width;
	//cout << score3 << endl;
	//cout << score << " " << score2 << endl;
	if (score2 == 0) score2 = 1;
	score = score * (double)(150 - (double)sedai) / 150 + score2 * sedai / 150;
	score = pow(score, 3) * score3;
	//cout << score * 1000 << endl;
	//cout << edge[0] << " " << edge[1] << endl;
	return score * 1000;
}




//�ӓ��m�̌�������////////////////////////////////////////////////////////////////////////////////////////////////////
bool CrossCheck(pieces p1, pieces p2) {
	double a, b, c, d;
	int ii, jj;
	for (int i = 0; i < p2.total_vertex; i++) {
		for (int j = 0; j < p1.total_vertex; j++) {
			if (i == p2.total_vertex - 1) ii = 0;
			else ii = i + 1;
			if (j == p1.total_vertex - 1) jj = 0;
			else jj = j + 1;
			a = (p2.vertex_x.at(i) - p2.vertex_x.at(ii)) * (p1.vertex_y.at(j) - p2.vertex_y.at(i)) + (p2.vertex_y.at(i) - p2.vertex_y.at(ii)) * (p2.vertex_x.at(i) - p1.vertex_x.at(j));
			b = (p2.vertex_x.at(i) - p2.vertex_x.at(ii)) * (p1.vertex_y.at(jj) - p2.vertex_y.at(i)) + (p2.vertex_y.at(i) - p2.vertex_y.at(ii)) * (p2.vertex_x.at(i) - p1.vertex_x.at(jj));
			c = (p1.vertex_x.at(j) - p1.vertex_x.at(jj)) * (p2.vertex_y.at(i) - p1.vertex_y.at(j)) + (p1.vertex_y.at(j) - p1.vertex_y.at(jj)) * (p1.vertex_x.at(j) - p2.vertex_x.at(i));
			d = (p1.vertex_x.at(j) - p1.vertex_x.at(jj)) * (p2.vertex_y.at(ii) - p1.vertex_y.at(j)) + (p1.vertex_y.at(j) - p1.vertex_y.at(jj)) * (p1.vertex_x.at(j) - p2.vertex_x.at(ii));
			if ((a * b <= 0) && (c * d <= 0)) return true;
		}
	}
	return false;
}


//���_�̓��O����////////////////////////////////////////////////////////////////////////////////////////////////////
bool InputCheck(pieces p1, pieces p2) {
	double deg = 0, inner_pro, count_deg = 0, pi = 3.141592653589793238;
	double vect[2][2];	//�x�N�g���𐶐� [i][j] => [i] �s��,��  [j] �x�N�g���̗v�f
	double abs_vect[2]; //�x�N�g���̐�Βl���i�[����ϐ�
	double vertexP2[2]; //�`�F�b�N���钸�_
	for (int j = 0; j < 1; j++) {
		vertexP2[0] = p2.vertex_x.at(j), vertexP2[1] = p2.vertex_y.at(j);

		for (int i = 0; i < p1.total_vertex; i++) {
			vect[0][0] = p1.vertex_x.at(i) - vertexP2[0], vect[0][1] = p1.vertex_y.at(i) - vertexP2[1]; //�x�N�g��1�𐶐�
			//�x�N�g��2�𐶐�
			if (i == p1.total_vertex - 1) {
				vect[1][0] = p1.vertex_x.at(0) - vertexP2[0];
				vect[1][1] = p1.vertex_y.at(0) - vertexP2[1];
			}
			else {
				vect[1][0] = p1.vertex_x.at(i + 1) - vertexP2[0];
				vect[1][1] = p1.vertex_y.at(i + 1) - vertexP2[1];
			}
			abs_vect[0] = sqrt((pow(vect[0][0], 2) + pow(vect[0][1], 2))); //�x�N�g��1�̑傫��
			abs_vect[1] = sqrt((pow(vect[1][0], 2) + pow(vect[1][1], 2))); //�x�N�g��2�̑傫��
			inner_pro = vect[0][0] * vect[1][0] + vect[0][1] * vect[1][1]; //�x�N�g��1,2�̓��ς��v�Z
			deg = inner_pro / (abs_vect[0] * abs_vect[1]); //�傫���Ɠ��ς���cos���Z�o
			//�O�ςɂ���]�������Z�o���C����������
			if (vect[0][0] * vect[1][1] - vect[0][1] * vect[1][0] <= 0) count_deg += acos(deg);
			else count_deg -= acos(deg);
			//1��]���Ă���Γ����ɑ���
			if (abs(count_deg) / (2 * pi) >= 0.9) {
				return true;
			}
		}
	}
	return false;
}

//�ˑR�ψ�////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<individual> Mutation(vector<individual> parent) {
	double mutation_per = 0.2;
	/*for (int i = 1; i < (int)parent.size(); i++) {
		if ((double)(rand() % 100) / 100 < mutation_per) parent.at(i).gen.at(2) = (double)(rand() % 365);
	}*/
	return parent;
}
