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
