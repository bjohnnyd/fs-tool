# TODO

1. Create all possible HLA structs at compile time.
  ```
     if you do it in build.rs, you don't need any serialization
     the data would not be stored in json or protobuf format but embedded into your tool
     in your build.rs, you would have code that
     1. reads the file and parses it (should be rather simple if it's tab-separated values)
     2. creates a bunch of rust code like const DATA : MyStruct { foo: 1234, bar: 5678, … }
     3. writes that to a file like generated/data.rs
  ```


