#include "diffexp2/json_codec.hpp"

#include "WolframLibrary.h"
#include "WolframNumericArrayLibrary.h"

#include <cstring>
#include <string>

namespace {

int return_bytes(WolframLibraryData lib_data, const std::string& value,
                 MArgument result) {
  const mint dims[1] = {static_cast<mint>(value.size())};
  MNumericArray bytes = nullptr;
  auto* api = lib_data->numericarrayLibraryFunctions;
  if (api == nullptr || api->MNumericArray_new(MNumericArray_Type_UBit8, 1, dims,
                                               &bytes) != LIBRARY_NO_ERROR) {
    return LIBRARY_MEMORY_ERROR;
  }
  if (!value.empty()) {
    std::memcpy(api->MNumericArray_getData(bytes), value.data(), value.size());
  }
  MArgument_setMNumericArray(result, bytes);
  return LIBRARY_NO_ERROR;
}

}  // namespace

extern "C" {

DLLEXPORT mint WolframLibrary_getVersion() { return WolframLibraryVersion; }
DLLEXPORT int WolframLibrary_initialize(WolframLibraryData) {
  return LIBRARY_NO_ERROR;
}
DLLEXPORT void WolframLibrary_uninitialize(WolframLibraryData) {}

DLLEXPORT int de2BackendInfo(WolframLibraryData lib_data, mint,
                             MArgument*, MArgument result) {
  return return_bytes(lib_data, diffexp2::backend_info_json(), result);
}

DLLEXPORT int de2RunRecurrence(WolframLibraryData lib_data, mint argc,
                               MArgument* args, MArgument result) {
  if (argc != 1) return LIBRARY_FUNCTION_ERROR;
  char* input = MArgument_getUTF8String(args[0]);
  std::string output;
  try {
    output = diffexp2::run_recurrence_json(input == nullptr ? "" : input);
  } catch (...) {
    if (input != nullptr) lib_data->UTF8String_disown(input);
    return LIBRARY_FUNCTION_ERROR;
  }
  if (input != nullptr) lib_data->UTF8String_disown(input);
  return return_bytes(lib_data, output, result);
}

}  // extern "C"
