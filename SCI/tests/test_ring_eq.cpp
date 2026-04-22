/*
 * Original Work copyright (c) 2021 Microsoft Research
 * Modified Work copyright (c) 2021 Microsoft Research
 *
 * Original Authors: Deevashwer Rathee, Mayank Rathee
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whome the Software is
 * furnished to do so, subject to the following conditions:
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
 * INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR
 * A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
 * OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 * Modified by Akash Shah
 */
#ifndef BATCHEQUALITY_H__
#define BATCHEQUALITY_H__
#include "OT/emp-ot.h"
#include "utils/emp-tool.h"
#include "Millionaire/bit-triple-generator.h"
#include <cmath>
#include <ctime>
#include <thread>
#include "BuildingBlocks/aux-protocols.h"
using namespace sci;
using namespace std;
#define WAN_EXEC
template <typename IO>
class BatchEquality
{
public:
	IO *io1 = nullptr;
	IO *io2 = nullptr;
	sci::OTPack *otpack1, *otpack2;
	TripleGenerator*triple_gen1, *triple_gen2;
	int party;
	int l, r, log_alpha, beta, beta_pow, batch_size, radixArrSize;
	int num_digits, num_triples_corr, num_triples_std, log_num_digits, num_cmps;
	int num_triples;
	uint8_t mask_beta, mask_r;
	Triple *triples_std;
	uint8_t *leaf_eq;
	uint8_t *digits;
	uint8_t **leaf_ot_messages;

	BatchEquality(int party,
				  int bitlength,
				  int log_radix_base,
				  int batch_size,
				  int num_cmps,
				  IOPack *iopack1,
				  IOPack *iopack2,
				  sci::OTPack *otpack1,
				  sci::OTPack *otpack2)
	{
		assert(log_radix_base <= 8);
		assert(bitlength <= 64);
		this->party = party;
		this->l = bitlength;
		this->beta = log_radix_base;
		this->io1 = iopack1->io;
		this->otpack1 = otpack1;
		this->io2 = iopack2->io;
		this->otpack2 = otpack2;
		this->batch_size = batch_size;
		this->num_cmps = num_cmps;
		this->triple_gen1 = new TripleGenerator(party, iopack1, otpack1);
		this->triple_gen2 = new TripleGenerator(3 - party, iopack2, otpack2);
		configure();
	}

	void configure()
	{
		this->num_digits = ceil((double)l / beta);
		this->r = l % beta;
		this->log_alpha = sci::bitlen(num_digits) - 1;
		this->log_num_digits = log_alpha + 1;
		this->num_triples = num_digits - 1;
		if (beta == 8)
			this->mask_beta = -1;
		else
			this->mask_beta = (1 << beta) - 1;
		this->mask_r = (1 << r) - 1;
		this->beta_pow = 1 << beta;
		this->triples_std = new Triple((num_triples)*batch_size * num_cmps, true);
	}

	~BatchEquality()
	{
		delete triple_gen1;
		delete triple_gen2;
	}

	void setLeafMessages(uint64_t *data)
	{

		if (this->party == sci::ALICE)
		{
			radixArrSize = batch_size * num_cmps;
		}
		else
		{
			radixArrSize = num_cmps;
		}
std::cout << "All Base OTs Done" << std::endl;
		digits = new uint8_t[num_digits * radixArrSize];
		leaf_eq = new uint8_t[num_digits * batch_size * num_cmps];

		for (int i = 0; i < num_digits; i++) // Stored from LSB to MSB
			for (int j = 0; j < radixArrSize; j++)
				if ((i == num_digits - 1) && (r != 0))
					digits[i * radixArrSize + j] = (uint8_t)(data[j] >> i * beta) & mask_r;
				else
					digits[i * radixArrSize + j] = (uint8_t)(data[j] >> i * beta) & mask_beta;

		if (party == sci::ALICE)
		{
			leaf_ot_messages = new uint8_t *[num_digits * num_cmps];
			for (int i = 0; i < num_digits * num_cmps; i++)
				leaf_ot_messages[i] = new uint8_t[beta_pow];

			// Set Leaf OT messages
			triple_gen1->prg->random_bool((bool *)leaf_eq, batch_size * num_digits * num_cmps);
			for (int i = 0; i < num_digits; i++)
			{
				for (int j = 0; j < num_cmps; j++)
				{
					if (i == (num_digits - 1) && (r > 0))
					{
#ifdef WAN_EXEC
						set_leaf_ot_messages(leaf_ot_messages[i * num_cmps + j], digits,
											 beta_pow, leaf_eq, i, j);
#else
						set_leaf_ot_messages(leaf_ot_messages[i * num_cmps + j], digits,
											 1 << r, leaf_eq, i, j);
#endif
					}
					else
					{
						set_leaf_ot_messages(leaf_ot_messages[i * num_cmps + j], digits,
											 beta_pow, leaf_eq, i, j);
					}
				}
			}
		}
	}

	void computeLeafOTs()
	{

		if (party == sci::ALICE)
		{

			// Perform Leaf OTs
#ifdef WAN_EXEC
			otpack1->kkot[beta-1]->send(leaf_ot_messages, num_cmps * (num_digits), 3);
#else
			if (r == 1)
			{
				otpack1->kkot[beta-1]->send(leaf_ot_messages, num_cmps * (num_digits - 1), 3);
				otpack1->iknp_straight->send(leaf_ot_messages + num_cmps * (num_digits - 1), num_cmps, 3);
			}
			else if (r != 0)
			{
				otpack1->kkot[beta-1]->send(leaf_ot_messages, num_cmps * (num_digits - 1), 3);
				if (r == 2)
				{
					otpack1->kkot[r-1]->send(leaf_ot_messages + num_cmps * (num_digits - 1), num_cmps, 3);
				}
				else if (r == 3)
				{
					otpack1->kkot[r-1]->send(leaf_ot_messages + num_cmps * (num_digits - 1), num_cmps, 3);
				}
				else if (r == 4)
				{
					otpack1->kkot[r-1]->send(leaf_ot_messages + num_cmps * (num_digits - 1), num_cmps, 3);
				}
				else
				{
					throw std::invalid_argument("Not yet implemented!");
				}
			}
			else
			{
				otpack1->kkot[0]->send(leaf_ot_messages, num_cmps * num_digits, 3);
			}
#endif
			// Cleanup
			for (int i = 0; i < num_digits * num_cmps; i++)
				delete[] leaf_ot_messages[i];
			delete[] leaf_ot_messages;
		}
		else // party = sci::BOB
		{	 // triple_gen1->generate(3-party, triples_std, _16KKOT_to_4OT);
		  //  Perform Leaf OTs
#ifdef WAN_EXEC
			otpack1->kkot[beta-1]->recv(leaf_eq, digits, num_cmps * (num_digits), 3);
#else
			if (r == 1)
			{
				otpack1->kkot[0]->recv(leaf_eq, digits, num_cmps * (num_digits - 1), 3);
				otpack1->iknp_straight->recv(leaf_eq + num_cmps * (num_digits - 1),
											 digits + num_cmps * (num_digits - 1), num_cmps, 3);
			}
			else if (r != 0)
			{
				otpack1->kkot[0]->recv(leaf_eq, digits, num_cmps * (num_digits - 1), 3);
				if (r == 2)
				{
					otpack1->kkot[r-1]->recv(leaf_eq + num_cmps * (num_digits - 1),
										  digits + num_cmps * (num_digits - 1), num_cmps, 3);
				}
				else if (r == 3)
				{
					otpack1->kkot[r-1]->recv(leaf_eq + num_cmps * (num_digits - 1),
										  digits + num_cmps * (num_digits - 1), num_cmps, 3);
				}
				else if (r == 4)
				{
					otpack1->kkot[r-1]->recv(leaf_eq + num_cmps * (num_digits - 1),
										   digits + num_cmps * (num_digits - 1), num_cmps, 3);
				}
				else
				{
					throw std::invalid_argument("Not yet implemented!");
				}
			}
			else
			{
				otpack1->kkot[0]->recv(leaf_eq, digits, num_cmps * (num_digits), 3);
			}
#endif

			// Extract equality result from leaf_res_cmp
			for (int i = 0; i < num_digits * num_cmps; i++)
			{
				for (int j = batch_size - 1; j >= 0; j--)
				{
					leaf_eq[j * num_digits * num_cmps + i] = (leaf_eq[i] >> j) & 1;
				}
			}
		}

		/*for(int i=0; i<10; i++) {
			for(int j=0;j<batch_size; j++) {
				std::cout<< (int)leaf_eq[j*num_digits*num_cmps+ i] << " ";
			}
			std::cout<< std::endl;
		}*/
		/*for (int i = 0; i < num_cmps; i++)
			res[i] = leaf_res_cmp[i];
 */
		// Cleanup
		delete[] digits;
	}

	void set_leaf_ot_messages(uint8_t *ot_messages,
							  uint8_t *digits,
							  int N,
							  uint8_t *mask_bytes,
							  int i,
							  int j)
	{
		for (int k = 0; k < N; k++)
		{
			ot_messages[k] = 0;
			for (int m = 0; m < batch_size; m++)
			{
				ot_messages[k] = ot_messages[k] | (((digits[i * radixArrSize + j * batch_size + m] == k) ^ mask_bytes[m * num_digits * num_cmps + i * num_cmps + j]) << m);
			}
		}
	}

	/**************************************************************************************************
	 *                         AND computation related functions
	 **************************************************************************************************/

	void generate_triples()
	{
		triple_gen2->generate(3 - party, triples_std, _16KKOT_to_4OT);
	}

	void traverse_and_compute_ANDs(uint8_t *res_shares)
	{
		// clock_gettime(CLOCK_MONOTONIC, &start);
		//  Combine leaf OT results in a bottom-up fashion
		int counter_std = 0, old_counter_std = 0;
		int counter_corr = 0, old_counter_corr = 0;
		int counter_combined = 0, old_counter_combined = 0;
		uint8_t *ei = new uint8_t[(num_triples * batch_size * num_cmps) / 8];
		uint8_t *fi = new uint8_t[(num_triples * batch_size * num_cmps) / 8];
		uint8_t *e = new uint8_t[(num_triples * batch_size * num_cmps) / 8];
		uint8_t *f = new uint8_t[(num_triples * batch_size * num_cmps) / 8];

		int old_triple_count = 0, triple_count = 0;

		for (int i = 1; i < num_digits; i *= 2)
		{
			int counter = 0;
			for (int j = 0; j < num_digits and j + i < num_digits; j += 2 * i)
			{
				for (int k = 0; k < batch_size; k++)
				{
					for (int m = 0; m < num_cmps; m += 8)
					{
						ei[(counter * batch_size * num_cmps + k * num_cmps + m) / 8] = triples_std->ai[(triple_count + counter * batch_size * num_cmps + k * num_cmps + m) / 8];
						fi[(counter * batch_size * num_cmps + k * num_cmps + m) / 8] = triples_std->bi[(triple_count + counter * batch_size * num_cmps + k * num_cmps + m) / 8];
						ei[(counter * batch_size * num_cmps + k * num_cmps + m) / 8] ^= sci::bool_to_uint8(leaf_eq + j * num_cmps + k * num_digits * num_cmps + m, 8);
						fi[(counter * batch_size * num_cmps + k * num_cmps + m) / 8] ^= sci::bool_to_uint8(leaf_eq + (j + i) * num_cmps + k * num_digits * num_cmps + m, 8);
					}
				}
				counter++;
			}
			triple_count += counter * batch_size * num_cmps;
			int comm_size = (counter * batch_size * num_cmps) / 8;

			if (party == sci::ALICE)
			{
				io1->send_data(ei, comm_size);
				io1->send_data(fi, comm_size);
				io1->recv_data(e, comm_size);
				io1->recv_data(f, comm_size);
			}
			else // party = sci::BOB
			{
				io1->recv_data(e, comm_size);
				io1->recv_data(f, comm_size);
				io1->send_data(ei, comm_size);
				io1->send_data(fi, comm_size);
			}

			for (int i = 0; i < comm_size; i++)
			{
				e[i] ^= ei[i];
				f[i] ^= fi[i];
			}

			counter = 0;
			for (int j = 0; j < num_digits and j + i < num_digits; j += 2 * i)
			{
				for (int k = 0; k < batch_size; k++)
				{
					for (int m = 0; m < num_cmps; m += 8)
					{
						uint8_t temp_z;
						if (party == sci::ALICE)
							temp_z = e[(counter * batch_size * num_cmps + k * num_cmps + m) / 8] & f[(counter * batch_size * num_cmps + k * num_cmps + m) / 8];
						else
							temp_z = 0;

						temp_z ^= f[(counter * batch_size * num_cmps + k * num_cmps + m) / 8] & triples_std->ai[(old_triple_count + counter * batch_size * num_cmps + k * num_cmps + m) / 8];
						temp_z ^= e[(counter * batch_size * num_cmps + k * num_cmps + m) / 8] & triples_std->bi[(old_triple_count + counter * batch_size * num_cmps + k * num_cmps + m) / 8];
						temp_z ^= triples_std->ci[(old_triple_count + counter * batch_size * num_cmps + k * num_cmps + m) / 8];
						sci::uint8_to_bool(leaf_eq + j * num_cmps + k * num_digits * num_cmps + m, temp_z, 8);
					}
				}
				counter++;
			}
			old_triple_count = triple_count;
		}

		for (int i = 0; i < num_cmps; i++)
		{
			res_shares[i] = 0;
			for (int j = 0; j < batch_size; j++)
			{
				res_shares[i] = res_shares[i] ^ leaf_eq[j * num_digits * num_cmps + i];
			}
		}

		// cleanup
		delete[] ei;
		delete[] fi;
		delete[] e;
		delete[] f;
	}
};

void computeLeafOTsThread(BatchEquality<NetIO> *compare)
{
	compare->computeLeafOTs();
}

void generate_triples_thread(BatchEquality<NetIO> *compare)
{
	compare->generate_triples();
}

void perform_batch_equality(uint64_t *inputs, BatchEquality<NetIO> *compare, uint8_t *res_shares)
{
	std::thread cmp_threads[2];
	compare->setLeafMessages(inputs);
	cmp_threads[0] = std::thread(computeLeafOTsThread, compare);
	cmp_threads[1] = std::thread(generate_triples_thread, compare);
	for (int i = 0; i < 2; ++i)
	{
		cmp_threads[i].join();
	}

	compare->traverse_and_compute_ANDs(res_shares);
}

#endif // BATCHEQUALITY_H__

int party, port = 32000;
string address = "127.0.0.1";

int dim = 1024;
int bw_x = 8;

sci::IOPack *iopackArr[4];
sci::OTPack *otpackArr[4];
AuxProtocols *aux;

int main(int argc, char **argv) {
  /************* Argument Parsing  ************/
  /********************************************/
  ArgMapping amap;
  amap.arg("r", party, "Role of party: ALICE = 1; BOB = 2");
  amap.arg("p", port, "Port Number");
  amap.arg("N", dim, "Number of ReLU operations");
  amap.arg("ip", address, "IP Address of server (ALICE)");

  amap.parse(argc, argv);


  /************ Generate Test Data ************/
  /********************************************/
  PRG128 prg;

  uint64_t *x = new uint64_t[dim];
  uint64_t *y = new uint64_t[dim];
  uint8_t *res = new uint8_t[dim];

  prg.random_data(x, dim * sizeof(uint64_t));

sci::IOPack *iopackArr[2];
sci::OTPack *otpackArr[2];
  /********** Setup IO and Base OTs ***********/
  /********************************************/
iopackArr[0] = new IOPack(party, port, address);
iopackArr[1] = new IOPack(party, port + 1, address);
otpackArr[0] = new OTPack(iopackArr[0], party);
otpackArr[1] = new OTPack(iopackArr[1], 3 - party);
aux = new AuxProtocols(party, iopackArr[0], otpackArr[0]);

  uint64_t total_comm = 0;
  uint64_t total_round = 0;
  total_comm = iopackArr[0]->get_comm();
  total_comm += iopackArr[1]->get_comm();
  total_round = iopackArr[0]->get_rounds();
  total_round += iopackArr[1]->get_rounds();

  auto start = clock_start();

  BatchEquality<NetIO>* compare;
  compare = new BatchEquality<NetIO>(party, bw_x, 2, 3, dim, iopackArr[0], iopackArr[1], otpackArr[0], otpackArr[1]);

  cout << "Bytes Sent\t" << iopackArr[0]->get_comm() + iopackArr[1]->get_comm() << " bytes" << endl;
  cout << "Rounds\t" << total_round << endl; 
  perform_batch_equality(x, compare, res);
  //aux->B2A(res, x, dim, 32);

  long long t = time_from(start);
  total_comm = iopackArr[0]->get_comm() - total_comm;
  total_comm += iopackArr[1]->get_comm();
  total_round = iopackArr[0]->get_rounds() - total_round;
  total_round += iopackArr[1]->get_rounds();

  cout << "Time\t" << t / (1000.0) << " ms" << endl;
  cout << "Bytes Sent\t" << total_comm << " bytes" << endl;
  cout << "Rounds\t" << total_round << endl;
  /******************* Cleanup ****************/
  /********************************************/
  delete[] x;
  delete[] y;
}