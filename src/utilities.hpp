/*-
 * SPDX-License-Identifier: BSD-2-Clause
 *
 * Copyright (c) 2021 Maarten L. Hekkelman, NKI-AVL
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <filesystem>
#include <regex>

#include <boost/program_options.hpp>

// --------------------------------------------------------------------

int a_main(int argc, char *const argv[]);
void print_what(const std::exception &ex);

// --------------------------------------------------------------------

boost::program_options::variables_map load_options(
	int argc, char *const argv[],
	boost::program_options::options_description &visible_options,
	boost::program_options::options_description &hidden_options,
	boost::program_options::positional_options_description &positional_options,
	const char *config_file_name);

// --------------------------------------------------------------------

class file_locator
{
  public:
	static void init(const std::filesystem::path &db_dir,
		const std::string &structure_name_pattern, const std::string &metadata_name_pattern)
	{
		s_instance.reset(new file_locator(db_dir, structure_name_pattern, metadata_name_pattern));
	}

	static std::filesystem::path get_structure_file(const std::string &id)
	{
		return s_instance->get_file(id, s_instance->m_structure_name_pattern);
	}

	static std::filesystem::path get_metdata_file(const std::string &id)
	{
		return s_instance->get_file(id, s_instance->m_metadata_name_pattern);
	}

  private:
	static std::unique_ptr<file_locator> s_instance;

	file_locator(const file_locator &) = delete;
	file_locator &operator=(const file_locator &) = delete;

	file_locator(const std::filesystem::path &db_dir,
		const std::string &structure_name_pattern, const std::string &metadata_name_pattern)
		: m_db_dir(db_dir)
		, m_structure_name_pattern(structure_name_pattern)
		, m_metadata_name_pattern(metadata_name_pattern)
	{
	}

	std::filesystem::path get_file(const std::string &id, std::string pattern)
	{
		std::string::size_type i;

		while ((i = pattern.find("${db-dir}")) != std::string::npos)
			pattern.replace(i, strlen("${db-dir}"), m_db_dir);

		std::regex rx(R"(\$\{id:(\d+):(\d+)\})");
		std::smatch m;

		while (std::regex_search(pattern, m, rx))
		{
			int start = stoi(m[1]);
			int length = stoi(m[2]);

			pattern.replace(m[0].first, m[0].second, id.substr(start, length));
		}

		while ((i = pattern.find("${id}")) != std::string::npos)
			pattern.replace(i, strlen("${id}"), id);

		return pattern;
	}

	const std::filesystem::path m_db_dir;
	const std::string m_structure_name_pattern;
	const std::string m_metadata_name_pattern;
};
