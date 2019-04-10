// Copyright 2019 Q1t BV
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

extern crate proc_macro;

use quote::quote;

macro_rules! gen_derive
{
    ($trait_name:ident, $fn_name:ident) => {
        #[proc_macro_derive($trait_name)]
        pub fn $fn_name(input: proc_macro::TokenStream) -> proc_macro::TokenStream
        {
            let ast: syn::DeriveInput = syn::parse(input).expect("Failed to build syntax tree");

            let name = &ast.ident;
            let gen = quote! {
                impl q1tsim::export::$trait_name for #name {}
            };
            gen.into()
        }
    }
}

gen_derive!(OpenQasm, open_qasm_derive);
gen_derive!(CQasm, c_qasm_derive);
gen_derive!(Latex, latex_derive);

#[proc_macro_derive(ExportGate)]
pub fn export_gate_derive(input: proc_macro::TokenStream) -> proc_macro::TokenStream
{
    let ast: syn::DeriveInput = syn::parse(input).expect("Failed to build syntax tree");

    let name = &ast.ident;
    let gen = quote! {
        impl q1tsim::export::OpenQasm for #name {}
        impl q1tsim::export::CQasm for #name {}
        impl q1tsim::export::Latex for #name {}
    };
    gen.into()
}
